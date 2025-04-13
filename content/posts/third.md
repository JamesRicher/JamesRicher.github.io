---
date: '2025-04-13'
draft: false
title: 'FFT-Based Water Rendering Part 3: Switching to the FFT'
params:
    math: true
categories: ["Water Rendering"]

cover:
    image: /images/fftWater_3/ocean_512.png
    alt: 'This is a post image'

ShowToc: true

weight: 3
---

{{< rawhtml >}}
<style>
    .pixelated {
        image-rendering: pixelated;
    }

    .smooth {
    image-rendering: smooth;
    }

    .pixelated {
    image-rendering: pixelated;
    }

    .crisp-edges {
    image-rendering: crisp-edges;
    }

    figcaption {
        font-style: italic;
        text-align: center;
    }
</style>
{{< /rawhtml >}}

At the moment, we have a fully functioning system for rendering our water. However, due to all of the additional DFTs we are performing to calculate our normals and displacement vectors, the simulation runs very slowly; using an \(N\) any higher than 256 results in single-digit frame rates on my machine. To remedy this and allow us to generate higher detail ocean surfaces we must replace the DFT with the much faster FFT. The particular FFT I will be implementing is the [Cooley–Tukey FFT](https://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm), for which we need to assume that \(N\) is always a power of two.

The full code for this project is available at my GitHub [here](https://github.com/JamesRicher/UnityFFTOcean).

## The 1D (Inverse) FFT
To develop an understanding of how the FFT works, we first consider the one-dimensional case. As a reminder, we can write the inverse DFT of a collection \(N\) frequency samples, \(F(0), \: F(1),\: ...\:, \:F(N - 1) \) as 

$$
f(n) = \sum_{k=0}^{N-1} f(k) \exp(2\pi i \frac{kn}{N})
$$

where \(n\) is an integer between \(0\) and \(N-1\) inclusive. This expression can be rewritten using the variable \(W_N\) called the "twiddle factor", defined as 

$$
W_N = \exp(\frac{2\pi i}{N})
$$

The twiddle factor has a few important properties we will need to apply later:
1. \(W_N^{k+N} = W_N^k\)
2. \(W_N^{k + \frac{N}{2}} = -W_N^k\)
3. \(W_N^{2k} = W_{\frac{N}{2}}^k\)

These properties are straightforward to derive from the definition \(W_N\).

Inserting \(W_N\) into the expression for the DFT, we find that

$$
f(n) = \sum_{k=0}^{N-1} f(k) W_N^{kn}
$$

Now comes the central idea of the FFT: we can split this DFT summation into two smaller DFTs. Each DFT will be of order \(\frac{N}{2}\) and run over either the odd indices or the even indices. This is done as follows:

$$
\begin{align}
    f(n) &= \sum_{k=0}^{N-1} f(k) W_N^{kn} \\
    &= \sum_{k=0}^{\frac{N}{2}-1} f(2k) W_N^{2kn} + \sum_{k=0}^{\frac{N}{2}-1} f(2k+1) W_N^{(2k + 1)n} \\
    &= \underbrace{\sum_{k=0}^{\frac{N}{2}-1} f(2k) W_{\frac{N}{2}}^{kn}}_{\frac{N}{2} \text{point DFT over even indices}} 
    + W_N^n \underbrace{\sum_{k=0}^{\frac{N}{2}-1} f(2k + 1) W_{\frac{N}{2}}^{kn}}_{\frac{N}{2} \text{point DFT over odd indices}} \\
\end{align}
$$

The above result is known as the **Danielson-Lanczos Lemma**. This exact same process can be repeated \(\log_2 (N) \) times for the smaller DFTs until the overall sum is no longer reducible. How this works in practice becomes much easier to see after studying some concrete cases for small \(N\). The FFT then builds on itself in a very natural and modular way for larger values of \(N\).

### The \(2\) Point FFT
The most basic case of the FFT is where \(N = 2\), which requires only one application of the Danielson-Lanczos Lemma. We can express this DFT as:

$$
\begin{align}
    f(n) &= \sum_{k=0}^1 F(k) W_2^{kn} \\
    &= \sum_{k=0}^0 F(2k) W_1^{kn} + W_2^n\sum_{k=0}^0 F(2k+1) W_1^{kn} \\
    &= F(0) + W_2^nF(1)
\end{align}
$$

Now, substituting concrete values for \(n\) into the above expression and using property \(2\) of \(W_N\) we find that:

$$
\begin{align}
    & f(0) = F(0) + W_2^0F(1) = F(0) + F(1) \\
    & f(1) = F(0) + W_2^1F(1) = F(0) - W_2^0F(1) = F(0) - F(1) 
\end{align}
$$

Sums of the form \(A + W_N^nB\) are called **butterfly operations** and are the fundamental building blocks of the FFT. These operations are neatly expressed in **butterfly diagrams**. For example, the butterfly diagram for the \(2\) point FFT is shown below.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_2.jpeg" alt="drawing" width="600"/>
            <figcaption>Butterfly diagram for \(2\) point FFT</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}


Butterfly diagrams can be understood as follows:
- Each value on the far left (\(F(0) \) and \(F(1)\)) represents an input to the system and enters through its own "channel". We can think of the inputs as flowing along their channel from left to right.
- A number, \(x\), above a channel indicates that the current contents of the channel should be multiplied by \(x\).
- Whenever a channel meets an outgoing arrow, its contents are added to the channel pointed to by that arrow. 
- The final outputs of the system are the values on the far right (\(f(0)\) and \(f(1)\)).

### The \(4\) Point FFT
Scaling by a power of two, we now consider the FFT for \(N=4\). Applying the Danielson-Lanczos Lemma once we obtain:

$$
\begin{align}
    f(n) &= \sum_{k=0}^3 F(k) W_4^{kn} \\
    &= \underbrace{\sum_{k=0}^1 F(2k) W_2^{kn}}_{N=2 \text{ DFT}} + W_4^n\underbrace{\sum_{k=0}^1 F(2k+1) W_2^{kn}}_{N=2 \text{ DFT}} \\
    &= E(n) + W_4^nO(n)
\end{align}
$$

Therefore, the first stage of the \(N=4\) FFT is to solve the two \(N=2\) DFTs, denoted \(E(n)\) for the even indices and \(O(n)\) for the odd indices. We have just seen exactly how to evaluate the FFT for \(N=2\) and so can apply that same logic here twice. Expressing \(f\) in terms of the subexpressions \(E\) and \(O\) and applying properties \(1\) and \(2\) of \(W_N\) we obtain:

$$
\begin{align}
    & f(0) = E(0) + W_4^0O(0)\\
    & f(1) = E(1) + W_4^1O(1)\\
    & f(2) = E(0) + W_4^2O(0) = E(0) - W_4^0O(0) \\
    & f(3) = E(1) + W_4^3O(1) = E(1) - W_4^1O(1) 
\end{align}
$$

The butterfly diagram for this FFT is shown below.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_4.jpeg" alt="drawing" width="700"/>
            <figcaption>Butterfly diagram for \(4\) point FFT</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_4_Condensed.jpeg" alt="drawing" width="700"/>
            <figcaption>Condensed butterfly diagram for \(4\) point FFT</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}



### The \(8\) Point FFT
Scaling by one final power of two, we now consider the FFT for \(N=8\) and proceed in exactly the same way.

$$
\begin{align}
    f(n) &= \sum_{k=0}^7 F(k) W_8^{kn} \\
    &= \underbrace{\sum_{k=0}^3 F(2k) W_4^{kn}}_{N=4 \text{ DFT}} + W_8^n\underbrace{\sum_{k=0}^3 F(2k+1) W_4^{kn}}_{N=4 \text{ DFT}} \\
    &= E(n) + W_8^nO(n)
\end{align}
$$

From the previous case, we have seen exactly how to efficiently solve the the two inner \(N=4\) DFTs. Again expressing \(f\) in terms of the subexpressions \(E\) and \(O\) (which now have a period of \(4\)) and applying properties \(1\) and \(2\) of \(W_N\) we obtain:

$$
\begin{align}
    & f(0) = E(0) + W_8^0O(0) \\
    & f(1) = E(1) + W_8^1O(1) \\
    & f(2) = E(2) + W_8^2O(2) \\
    & f(3) = E(3) + W_8^3O(3) \\
    & f(4) = E(0) + W_8^4O(0) = E(0) - W_8^0O(0)\\
    & f(5) = E(1) + W_8^5O(1) = E(1) - W_8^1O(1) \\
    & f(6) = E(2) + W_8^6O(2) = E(2) - W_8^2O(2) \\
    & f(7) = E(3) + W_8^7O(3) = E(3) - W_8^3O(3)
\end{align}
$$

The butterfly diagram for this FFT is shown below.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_8.jpeg" alt="drawing" width="700"/>
            <figcaption>Butterfly diagram for \(8\) point FFT</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_8_Condensed.jpeg" alt="drawing" width="700"/>
            <figcaption>Condensed butterfly diagram for \(8\) point FFT</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

### Bit Reversal
It is is important to note that the order of the inputs to the FFT is not the original order. This is due to the way the algorithm recursively breaks down the inputs into odd and even components, effectively splitting the input samples into pairs sharing the same congruency class modulo \(\frac{N}{2}\). The correct ordering can be neatly found by bit-reversing the index of each input (i.e. mirroring its binary representation) and using this reversed value as its new input index in the FFT. For example, the following is the bit-reversal table for \(N=8\):

index | binary rep. | reversed binary rep.| new index 
--- | --- | --- | ----
0 | 000 | 000 | 0
1 | 001 | 100 | 4
2 | 010 | 010 | 2
3 | 011 | 110 | 6
4 | 100 | 001 | 1
5 | 101 | 101 | 5
6 | 110 | 011 | 3
7 | 111 | 111 | 7


## Implementing the 1D FFT
Putting this all together, the FFT can be implemented in C# as follows:

``` cs
private Complex[] IFFT(Complex[] samples)
{
    int N = samples.Length;
    int stages = (int)Mathf.Log(N,2);
    Complex[] output = BitReverseArray(samples, stages);

    // For each stage
    for (int i = 1; i <= stages; i++)
    {
        int m = (int)Mathf.Pow(2f, i);

        Complex omegaM = Complex.CreateFromExponential(1f, 2*Mathf.PI/m);

        // For each group in the current stage
        for (int k = 0; k < N; k += m)
        {
            // For each pair in the group
            Complex omega = Complex.CreateFromParts(1f,0f);
            for (int j = 0; j < m/2; j++)
            {
                Complex t = omega * output[k+j + (m/2)];
                Complex u = output[k+j];

                output[k+j] = u + t;
                output[k+j + (m/2)] = u - t;
                omega *= omegaM;
            }
        }
    }

    return output.Select(x => x/(float)N).ToArray();
}
```

## The 2D FFT
Now, moving back towards the original problem, we consider how we can apply the FFT to a 2D DFT. Taking the spatial domain to be \(\mathbf{x} = (x, y)\) and the frequency domain to be \(\mathbf{k} = (k, j)\), we can write and rearrange a standard 2D (inverse) DFT as follows:

$$
\begin{align}
    f(\mathbf{x}) &= \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N} \mathbf{k} \cdot \mathbf{x}) \\
    &= \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(k,j) \exp (\frac{2\pi i}{N} (kx + jy)) \\
    &= \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(k,j) \exp (\frac{2\pi i}{N}kx) \exp (\frac{2\pi i}{N}jy) \\
    &= \sum_{j=0}^{N-1} \underbrace{\Biggl (\sum_{k=0}^{N-1} f(k,j) \exp (\frac{2\pi i}{N}kx) \biggr)}_{N \text{ point 1D DFT}}  \exp (\frac{2\pi i}{N}jy) 
\end{align}
$$

As we can see, for each constant value of \(j\) the underlined term represents a standard 1D DFT in the horizontal direction which can be solved using the FFT. The outer expression can then be seen to be another series of \(N\) 1D DFTs, but this time in the vertical direction. This is the approach we will take to turn our spectrum texture into a heightfield, as well as to compute the normals and displacement vectors. We will proceed in three stages:

1. Compute \(N\) FFTs, taking each row of \(N\) points in the spectrum texture as the inputs. Write this output to a texture.
2. Compute \(N\) more FFTs, taking each column of \(N\) points in the output texture of step \(1\) as the inputs.
3. Multiply each output pixel by either \(1\) or \(-1\) to account for the fact that our frequency domain ranges over \(-\frac{N}{2}\) to \(\frac{N}{2} - 1\) rather than \(0\) to \(N-1\). Also, scale the output by \(\frac{1}{N^2}\) as is standard.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_Step1.png" alt="drawing" width="500"/>
            <figcaption>2D DFT step 1</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/FFT_Step2.png" alt="drawing" width="500"/>
            <figcaption>2D DFT step 2</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

## The FFT Using Compute Shaders

The final hurdle to overcome in this section is how we actually implement the FFT on the GPU. This is more complicated than the CPU implementation presented earlier as, on the GPU, we must compute each output pixel in isolation. My implementation of the FFT on the GPU is an based on the implementation by Fynn-Jorin Flügge in [Realtime GPGPU FFT Ocean Water Simulation](https://www.google.com/search?client=safari&rls=en&q=Realtime+GPGPU+FFT+Ocean+Water+Simulation&ie=UTF-8&oe=UTF-8&safe=active).

### The Butterfly Texture
The key ingredient to the GPU FFT is the **butterfly texture**. For an \(N\) point FFT, the butterfly texture is a texture of width \(\log_2(N)\) (i.e. the number of stages in the FFT) and height \(N\). Each pixel in this texture can be thought of as one of the \(N\log_2(N)\) butterfly operations in the FFT, with the horizontal position representing the stage and the vertical position representing the channel. The contents of each pixel contain the information required to compute the corresponding butterfly operation. This information is stored in the texture's colour channels as follows:

- red - real part of the relevant twiddle factor.
- green - imaginary part of the relevant twiddle factor. 
- blue - index of the first sample.
- alpha - index of the second sample (the one to be multiplied by the twiddle factor).

As an example, the following table represents the contents of the butterfly texture for \(N=8\), with each cell corresponding to a pixel as an (R, G, B, A) vector:

{{< rawhtml >}}
    <table>
        <tr>
            <td>(\(-W_2^0, 3, 7\))</td>
            <td>(\(-W_4^1, 5, 7\))</td>
            <td>(\(-W_8^3, 3, 7\))</td>
        </tr>
        <tr>
            <td>(\(W_2^0, 3, 7\))</td>
            <td>(\(-W_4^0, 4, 6\))</td>
            <td>(\(-W_8^2, 2, 6\))</td>
        </tr>
        <tr>
            <td>(\(-W_2^0, 1, 5\))</td>
            <td>(\(W_4^1, 5, 7\))</td>
            <td>(\(-W_8^1, 1, 5\))</td>
        </tr>
        <tr>
            <td>(\(W_2^0, 1, 5\))</td>
            <td>(\(W_4^0, 4, 6\))</td>
            <td>(\(-W_8^0, 0, 4\))</td>
        </tr>
        <tr>
            <td>(\(-W_2^0, 2, 6\))</td>
            <td>(\(-W_4^1, 1, 3\))</td>
            <td>(\(W_8^3, 3, 7\))</td>
        </tr>
        <tr>
            <td>(\(W_2^0, 2, 6\))</td>
            <td>(\(-W_4^0, 0, 2\))</td>
            <td>(\(W_8^2, 2, 6\))</td>
        </tr>
        <tr>
            <td>(\(-W_2^0, 0, 4\))</td>
            <td>(\(W_4^1, 1, 3\))</td>
            <td>(\(W_8^1, 1, 5\))</td>
        </tr>
        <tr>
            <td>(\(W_2^0, 0, 4\))</td>
            <td>(\(W_4^0, 0, 2\))</td>
            <td>(\(W_8^0, 0, 4\))</td>
        </tr>
    </table>
{{< /rawhtml >}}

Here are a couple of actual butterfly textures. It is important to note that the vertical orientation of these textures is flipped when compared to standard butterfly diagrams due to texture coordinates starting from the bottom left.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_3/butterfly_32.png" alt="drawing" width="90"/>
            <figcaption>\(N = 32\) butterfly texture shown at its actual resolution</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_3/butterfly_128_stretched.png" alt="drawing" width="500"/>
            <figcaption>\(N=128\) butterfly texture stretched to be more readable</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

The following is a Unity compute shader that will generate the butterfly texture for a given choice of \(N\):

```hlsl
#pragma kernel GenerateButterflyTexture

#include "Assets/ShaderIncludes/ComplexNumbers.cginc"

RWTexture2D<float4> _butterflyTexture;
StructuredBuffer<int> _bitReversedIndices;

float _PI; // globally defined shader constant
int _N;

[numthreads(1,8,1)]
void GenerateButterflyTexture (uint3 id : SV_DispatchThreadID)
{
    // id.x = current stage (from zero)
    // id.y = current input index
    // ouptut (x,y) = twiddle factor, (z,w) = input indices for butterfly operation
    // Note: this texture is upside down when compared to conventional butterfly diagrams

    // Getting twiddle factors
    int currentStage = id.x + 1;
    int groupCount = _N >> currentStage;
    float k = id.y * groupCount % _N;

    float ang = k * 2 * _PI / (float)_N;
    float2 W = ComplexExp(ang);

    // Getting indices 
    uint butterflySpan = (uint)pow(2, id.x);
    uint groupSize = (uint)pow(2, id.x + 1);
    bool inTopButterfly = ((id.y % groupSize) < butterflySpan);
    bool inStageOne = id.x == 0;

    int2 indices = inTopButterfly ? int2(id.y, id.y + butterflySpan) : int2(id.y - butterflySpan, id.y);
    indices = lerp(indices, int2(_bitReversedIndices[indices.x], _bitReversedIndices[indices.y]), inStageOne);

    _butterflyTexture[id.xy] = float4(W.x, W.y, indices.x, indices.y);
}

```

### The Butterfly Operations
Using the butterfly texture, a horizontal butterfly operation can be computed as follows:

```hlsl
float4 butterflyData = _butterflyTexture[int2(_currentStage, id.x)];
float2 W = butterflyData.xy;

float2 p = _inputTexture[int2(butterflyData.z, id.y)].xy;
float2 q = _inputTexture[int2(butterflyData.w, id.y)].xy;

float2 output = p + ComplexMult(W,q);
```

where `_currentStage` is an `int` representing the current FFT stage and `id` is an `int2` representing the current pixel ID being worked on. 

As each input element is used in two distinct butterfly calculations, we cannot overwrite it in the input texture with `output`. Therefore, we must write the output of a stage of the FFT to its own texture. In practice, we use two textures called `_pingpong0Texture` and `_pingpong1Texture`, one of which acts as the input and the other as the output, depending on the value of an integer `_pingpong`. When `_pingpong == 0`, `_pingpong0Texture` is the input and when `_pingpong == 1`, `_pingpong1Texture` is the input. The value of this parameter is flipped on the CPU after each FFT step and passed to the compute shader before it is dispatched. 

Here is a complete Unity compute shader for a horizontal butterfly operation using the pingpong textures:

```hlsl
#pragma kernel HorizontalButterfly

#include "Assets/ShaderIncludes/ComplexNumbers.cginc"

Texture2D<float4> _butterflyTexture;
RWTexture2D<float4> _pingpong0Texture;
RWTexture2D<float4> _pingpong1Texture;

int _N;
float _PI;
int _currentStage;
int _pingpong; 

[numthreads(8,8,1)]
void HorizontalButterfly(uint3 id : SV_DispatchThreadID)
{
    float4 butterflyData = _butterflyTexture[int2(_currentStage, id.x)];
    float2 W = butterflyData.xy;

    if (_pingpong == 0)
    {
        float2 p = _pingpong0Texture[int2(butterflyData.z, id.y)].xy;
        float2 q = _pingpong0Texture[int2(butterflyData.w, id.y)].xy;

        float2 output = p + ComplexMult(W,q);
        _pingpong1Texture[id.xy] = float4(output.x, output.y, 0, 0);
    }
    else
    {
        float2 p = _pingpong1Texture[int2(butterflyData.z, id.y)].xy;
        float2 q = _pingpong1Texture[int2(butterflyData.w, id.y)].xy;

        float2 output = p + ComplexMult(W,q);
        _pingpong0Texture[id.xy] = float4(output.x, output.y, 0, 0);
    }
}
```

### Permuting and Scaling
The final step in the FFT is to permute the output to reflect the fact that our frequency domain centred about zero, rather than starting from zero. What we are calculating at the moment is 

$$
\begin{align}
    & \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N} \mathbf{k} \cdot \mathbf{x}) \\
    &= \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N}kx)\exp (\frac{2\pi i}{N}jy)
\end{align}
$$

Rearranging this into the desired form we obtain:

$$
\begin{align}
    & \sum_{k=0}^{N-1} \sum_{j=0}^{N-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N}kx)\exp (\frac{2\pi i}{N}jy)\\
    &= \sum_{k=-\frac{N}{2}}^{\frac{N}{2}-1} \sum_{j=-\frac{N}{2}}^{\frac{N}{2}-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N}(k+ \frac{N}{2})x)\exp (\frac{2\pi i}{N}(j + \frac{N}{2})y) \\
    &= \exp(\pi i)^x \exp(\pi i)^y  \sum_{k=-\frac{N}{2}}^{\frac{N}{2}-1} \sum_{j=-\frac{N}{2}}^{\frac{N}{2}-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N}kx)\exp (\frac{2\pi i}{N}jy) \\
    &= (-1)^{x+y}  \sum_{k=-\frac{N}{2}}^{\frac{N}{2}-1} \sum_{j=-\frac{N}{2}}^{\frac{N}{2}-1} f(\mathbf{k}) \exp (\frac{2\pi i}{N}kx)\exp (\frac{2\pi i}{N}jy) \\
\end{align}
$$

(Note that the argument \(\mathbf{k}\) of \(f\) does not change in the above expression as this argument is purely symbolic). Therefore, to obtain the final heightfield, we must multiply each pixel output by \((-1)^{x+y}\) where \(x\) and \(y\) are the pixel's ID. At this stage we also scale the output by the standard normalisation factor \(\frac{1}{N^2}\).

### Putting it All Together
With all of the elements of our FFT written, we can now put them together in C#. We dispatch the horizontal butterfly compute shader \(\log_2(N)\) times to complete the \(N\) horizontal FFTs, setting the stage and pingpong variables before each dispatch. We then do the same in the vertical direction using the vertical butterfly compute shader. The output of this is then permuted and scaled to create the final output texture. 

This can be implemented in C# as follows:

```cs
public void IFFT2D(RenderTexture input, RenderTexture output)
{
    _pingpong = 0;
    Graphics.Blit(input, _pingpong0RT);

    for (int i = 0; i < _totalStages; i++)
    {
        _fftCompute.SetInt("_currentStage", i);
        _fftCompute.SetInt("_pingpong", _pingpong);
        _fftCompute.Dispatch(_horButterflyKernelID, _workGroupCount, _workGroupCount, 1);
        _pingpong = (_pingpong + 1) % 2;
    }

    for (int i = 0; i < _totalStages; i++)
    {
        _fftCompute.SetInt("_currentStage", i);
        _fftCompute.SetInt("_pingpong", _pingpong);
        _fftCompute.Dispatch(_verButterflyKernelID, _workGroupCount, _workGroupCount, 1);
        _pingpong = (_pingpong + 1) % 2;
    }

    _permuteAndScaleCompute.SetTexture(_permuteAndScaleKernelID, "_outputTexture", output);
    _permuteAndScaleCompute.Dispatch(_permuteAndScaleKernelID, _workGroupCount, _workGroupCount, 1);
}
```

## Rendering With Higher \(N\)
Making the jump from the DFT to the FFT, we can now render the water's surface at higher resolutions, being able to push \(N\) up to \(1024\) on my machine. The original DFT implementation required approximately \(N^2\) operations to compute the heightfield, which for \(N= 1024\) means \(1048576\) total operations. The FFT on the other hand requires approximately \(N\log_2(N)\) operations, which for the same value of \(N\) means only \(10240\) total operations. As a result, for \(N=1024\) the FFT is over \(100\) times faster than the DFT.
{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_3/ocean_1024.jpg" alt="drawing" width="800"/>
            <figcaption>Ocean surface using \(N=1024\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

## References
- [Fynn-Jorin Flügge: Realtime GPGPU FFT Ocean Water Simulation](https://www.google.com/search?client=safari&rls=en&q=Realtime+GPGPU+FFT+Ocean+Water+Simulation&ie=UTF-8&oe=UTF-8&safe=active)
- [Always Learn: A DFT and FFT Tutorial](https://www.alwayslearn.com/DFT%20and%20FFT%20Tutorial/DFTandFFT_FFT_Overview.html)
- Ronald N. Bracewell: The Fourier Transform and its Applications (book)