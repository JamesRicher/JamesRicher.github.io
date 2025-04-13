---
date: '2025-04-08'
draft: false
title: 'FFT-Based Water Rendering Part 2: Normals and Ocean Swell'
url: 'FFTWaterPart2'
params:
    math: true
showToc: true

tags: ["Unity", "HLSL", "C#", "Maths"]

cover:
    image: /images/fftWater_2/swell2.png
    alt: 'This is a post image'

categories: ["Water Rendering"]

weight: 2
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

Following on from the first part in this series, we now have our ocean heightfield but no way to make out much of the wave details as we are not currently computing the surface normals. In this part, I will show how these normals can be calculated, as well as improving the waves' presentation by adding some vertex movement in the \(xz\) plane to simulate ocean swell. 

The full code for this project is available at my GitHub [here](https://github.com/JamesRicher/UnityFFTOcean).

## Computing the Normals
There are several ways to go about getting the surface normals. One way is to approximate the normals using a central difference method to find estimates for the partial derivatives of the surface. While this is fast, we can actually do better and obtain an exact expression for the normal at each point using the linearity of the DFT. Treating our ocean as a surface, \(f\), parameterised by \(\mathbf{x} = (x,z)\) we have that

$$
f(x,z) = \begin{bmatrix}  
            x \\
            h(x,z) \\
            z
            \end{bmatrix}
$$
Taking the partial derivatives to find the vectors spanning to tangent plane and taking their cross product, we compute the normal, \(N(\mathbf{x}) \), as

$$
N(\mathbf{x}) = \frac{\partial f}{\partial x} \times \frac{\partial f}{\partial z} =
            \begin{bmatrix}  
            1 \\
            \frac{\partial h}{\partial x} \\
            0
            \end{bmatrix} \times
            \begin{bmatrix}  
            0 \\
            \frac{\partial h}{\partial z} \\
            1
            \end{bmatrix} =
            \begin{bmatrix}  
            \frac{\partial h}{\partial x} \\
            -1 \\
            \frac{\partial h}{\partial z}
            \end{bmatrix}
$$
Note that this normal will always have a negative \(y\) component, therefore the normal we are actually looking for is \(-N(\mathbf{x})\) (and of course scaled to be of unit length, \(\frac{-N(\mathbf{x})}{||N(\mathbf{x})||}\)).

Computing these partial derivatives using the linearity of the Fourier Transform is relatively straightforward, however we must be aware of some subtleties. It turns out that we cannot simply take our partial derivatives to be

$$
\begin{align}
\frac{\partial h}{\partial x} = \sum_{\mathbf{k}} ik_x\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) \\
\frac{\partial h}{\partial z} = \sum_{\mathbf{k}} ik_z\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) \\
\end{align}
$$

and must instead use

$$
\begin{align}
\frac{\partial h}{\partial x} = \sum_{\mathbf{k}} d_x(k_x)\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) \\
\frac{\partial h}{\partial z} = \sum_{\mathbf{k}} d_z(k_z)\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) \\
\end{align}
$$

where

$$
\begin{align}
d_x(k_x) = \begin{cases} 
            0, & \text{if \(k_x = -\frac{N\pi}{L_x}\)}\\
            ik_x, & \text{otherwise}
        \end{cases} \\
d_z(k_z) = \begin{cases} 
            0, & \text{if \(k_z = -\frac{M\pi}{L_z}\)}\\
            ik_z, & \text{otherwise}
        \end{cases}
\end{align}
$$

The reasoning behind this is discussed by Steven Johnson in [\'Notes on FFT-based differentiation\'](https://math.mit.edu/%7Estevenj/fft-deriv.pdf).

In implementation, we can write the functions \(d_x(k_x)\tilde{h}(\mathbf{k},t)\) and \(d_z(k_z)\tilde{h}(\mathbf{k},t)\) to two textures and apply a regular (inverse) DFT to each to obtain the partial derivatives. 

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_2/xslope.png" alt="drawing" width="500"/>
            <figcaption>\( \frac{\partial h}{\partial x}\) packed into the red channel</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_2/zslope.png" alt="drawing" width="500"/>
            <figcaption>\( \frac{\partial h}{\partial z}\) packed into the red channel</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

## Optimising the Fourier Transforms
Note that, as with the original spectrum \(\tilde{h}(\mathbf{k}, t)\), the two functions we take the (inverse) DFT of to obtain \(\frac{\partial h}{\partial x}\) and \(\frac{\partial h}{\partial x}\) are also hermitian and therefore yield purely real-valued functions. We can use this property to our advantage in order to half the number of (inverse) DFTs we are performing. 

Consider two arbitrary discrete hermitian spectra \(\tilde{f}(\mathbf{k})\) and \(\tilde{g}(\mathbf{k})\). We have then that

$$
\begin{align}
f(\mathbf{x}) = \mathcal{F}^{-1}[\tilde{f}(\mathbf{k})] \in \mathbb{R} \\
g(\mathbf{x}) = \mathcal{F}^{-1}[\tilde{g}(\mathbf{k})] \in \mathbb{R} 
\end{align}
$$

Now consider the spectrum \(\tilde{c}(\mathbf{k})\) obtained by combining \(\tilde{f}\) and \(\tilde{g}\) as follows:

$$
\tilde{c}(\mathbf{k}) = \tilde{f}(\mathbf{k}) + i \tilde{g}(\mathbf{k})
$$

Then the IFT of this spectrum can be found as 

$$  
\begin{align}
    c(\mathbf{x}) &= \mathcal{F}^{-1}[\tilde{c}(\mathbf{k})] \\
    &= \sum_{\mathbf{k}} (\tilde{f}(\mathbf{k}) + i \tilde{g}(\mathbf{k}))\exp (i \mathbf{k} \cdot \mathbf{x}) \\
    &= \sum_{\mathbf{k}}\tilde{f}(\mathbf{k})\exp (i \mathbf{k} \cdot \mathbf{x}) + i\sum_{\mathbf{k}}\tilde{g}(\mathbf{k})\exp (i \mathbf{k} \cdot \mathbf{x}) \\
    &= f(\mathbf{x}) + ig(\mathbf{x})
\end{align}
$$

The upshot of this is that we can calculate \(f(\mathbf{x})\) and \(g(\mathbf{x})\) using only **one** (inverse) DFT using the fact that \(f(\mathbf{x}) =\operatorname{Re}(c(\mathbf{x}))\) and \(g(\mathbf{x}) =\operatorname{Im}(c(\mathbf{x})) \). Therefore, anytime we have a pair of hermitian spectra, we can pack the corresponding spectra in one texture as above, apply a single (inverse) DFT, and then extract the individual output functions from the red and green channels of the resulting texture.

## Adding Swell
One of the most glaringly unrealistic elements of our ocean surface at the moment is how rounded our wave peaks appear. In a realistic ocean simulation, we would expect to see a larger amount of variation, with much sharper wave peaks and broader troughs, creating the kind of "swell" exhibited by a real ocean surface. To achieve this, we can additionally displace our vertices in the \(xz\) plane, moving them towards the steeper areas of the heightfield. Practically, this can be done by displacing the vertices using the 2D displacement vector \(\mathbf{D}(\mathbf{x},t) = (D_x(\mathbf{x},t), D_z(\mathbf{x},t))\) defined by

$$
\begin{align}
    D_x(\mathbf{x},t) = \sum_{\mathbf{k}} d_x(k_x, k)\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) \\
    D_z(\mathbf{x},t) = \sum_{\mathbf{k}} d_z(k_z, k)\tilde{h}(\mathbf{k},t)\exp (i\mathbf{k} \cdot \mathbf{x}) 
\end{align}
$$

where

$$
\begin{align}
d_x(k_x, k) = \begin{cases} 
            0, & \text{if \(k_x = -\frac{N\pi}{L_x}\) or \(k =0\)}\\
            \frac{ik_x}{k}, & \text{otherwise}
        \end{cases} \\
d_z(k_z, k) = \begin{cases} 
            0, & \text{if \(k_z = -\frac{M\pi}{L_z}\) or \(k =0\)}\\
            \frac{ik_z}{k}, & \text{otherwise}
        \end{cases}
\end{align}
$$

Roughly speaking, we can think of \(\mathbf{D}(\mathbf{x},t)\) as the gradient vector of the heightfield and therefore pointing in the direction of its steepest slope. 

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_2/xdisp.png" alt="drawing" width="500"/>
            <figcaption>\( D_x\) packed into the red channel</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_2/zdisp.png" alt="drawing" width="500"/>
            <figcaption>\( D_z\) packed into the red channel</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}
Combining this with the heightfield, we can express the final displaced position of the vertex at \(\mathbf{x}\) as 

$$
\begin{bmatrix}
    x + \lambda D_x(\mathbf{x}, t) \\
    h(\mathbf{x}, t) \\
    z + \lambda D_z(\mathbf{x}, t)
\end{bmatrix}
$$

where \(\lambda \in \mathbb{R}, \lambda \geq 0\) is a parameter that controls the strength of the swell. \(\lambda\) should be kept relatively small to avoid "loops" forming, where vertices are displaced past the wave crest in each direction.

## Computing the Modified Normals
As we have modified the surface, we need to recalculate the normals to reflect this. The new correct expression for the normal is given by

$$
\begin{bmatrix}  
    -\frac{\frac{\partial h}{\partial x}}{1 + \lambda \frac{\partial D_x}{\partial x}} \\
    1 \\
    -\frac{\frac{\partial h}{\partial z}}{1 + \lambda \frac{\partial D_z}{\partial z}}
\end{bmatrix}
$$

(scaled to unit length). The partial derivatives \(\frac{\partial D_x}{\partial x}\) and \(\frac{\partial D_z}{\partial z}\) needed for the above expression can be computed in exactly the same way as we found \( \frac{\partial h}{\partial x} \) and \(\frac{\partial h}{\partial z}\).

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_2/no-swell.png" alt="drawing" width="1000"/>
            <figcaption>Ocean surface at \(N = 256\) and \(\lambda = 0\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

{{< rawhtml >}}
    <center>
        <figure>
            <img class="smooth" src="/images/fftWater_2/swell1.png" alt="drawing" width="1000"/>
            <figcaption>Ocean surface at \(N = 256\) and \(\lambda = 1\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

For now, we just render the ocean surface using a very basic shader that only takes into account skybox reflections modulated by a Fresnel factor, however we will come back and improve this shading model later.


## References
- [Thomas Gamper: Ocean Surface Generation and Rendering](https://www.cg.tuwien.ac.at/research/publications/2018/GAMPER-2018-OSG/GAMPER-2018-OSG-thesis.pdf), from which I took the expressions for \(\mathbf{D}(\mathbf{x}, t)\) and normals taking into account swell.
- [Steven Johnson: Notes on FFT-based differentiation](https://math.mit.edu/%7Estevenj/fft-deriv.pdf).
- [Jerry Tessendorf: Simulating Ocean Water](https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf)