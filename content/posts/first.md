---
date: '2025-04-04'
draft: false
title: 'FFT-Based Water Rendering Part 1: Building the Heightfield'
tags: ["Unity", "HLSL", "Maths"]
categories: ["Water Rendering"]

cover:
    image: /images/fftWater_1/part1_cover.png
    alt: 'This is a post image'

params:
    math: true
showToc: true

summary: write a summary here

weight: 1
---

{{< rawhtml >}}
<style>
    .pixelated {
        image-rendering: pixelated;
    }

    figcaption {
        font-style: italic;
        text-align: center;
    }
</style>
{{< /rawhtml >}}


**Welcome to my first blog post!** In this short series, I will be walking through a method of ocean surface rendering centred around the FFT, with an implementation in Unity. This approach was first developed by Jerry Tessendorf in the 2001 paper [\'Simulating Ocean Water\'](https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf) and has since been used in films such as Waterworld and Titanic, and more recently for real-time applications such as in Rare's Sea Of Thieves.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/sea-of-thieves.png" alt="drawing" width="800"/>
            <figcaption>Sea of Thieves</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

The full code for this project is available at my GitHub [here](https://github.com/JamesRicher/UnityFFTOcean).

## Tessendof's Method
The general approach of this technique is the same as most more lightweight methods of water rendering - sum a large number of sinusoids per-vertex with varying speed, direction and wavelength until you have a convincingly unpredictable heightfield. [GPU Gems](https://developer.nvidia.com/gpugems/gpugems/part-i-natural-effects/chapter-1-effective-water-simulation-physical-models) has a great introductory article on this which goes into more detail.

The primary difference with this method however is that we do all of our calculations directly in the frequency domain, and "pull out" the heightfield using the Discrete Fourier Transform (DFT). There are two key reasons that this approach is so powerful:
- By leveraging the power of the Fast Fourier Transform (FFT), we can sum many more waves at interactive speeds than if we were working in the spatial domain.
- It is particularly easy to apply real world ocean data when working in the frequency domain by making use of oceanographic spectra.

There are two key pairs of parameters we use to control the system. Firstly, we can change \( L_x \) and \(L_z\) which represent the physical size of the patch of water we are rendering (in metres) in each dimension. Secondly, we can change the integer parameters \(N\) and \(M\). These integers represent the resolution of the frequency space in each dimension and, by extension, the resolution of the generated heightfield. For my implementation I always take \(N = M = 2^r \) for some \(r \in \mathbb{N}\) to make the FFT implementation more straightforward, with values ranging from \(128\) to \(1024\). As a rough guide, Tessendorf states that the physical resolutions \( dx = \frac{L_x}{N}\) and \( dz = \frac{L_z}{M}\) need not go below \(2 \mathrm{cm} \).

## The Phillip's Spectrum
For my implementation, I have chosen to use the Phillip's Spectrum, \( P(\mathbf{k})\), as my oceanographic spectrum. This function takes the form 

$$
P(\mathbf{k}) = \frac{A * \exp(\frac{-1}{(kL)^2})) | \hat{\mathbf{k}} \cdot \hat{\mathbf{\omega}} |^2}{k^4}
$$
where:
- \( A \in \mathbb{R}\) is a multiplicative constant.
- \(\mathbf{k} \) is a wave vector with magnituded \(k\) and \(\hat{\mathbf{k}} = \frac{\mathbf{k}}{k} \).
- \(L = \frac{V^2}{g}\) where \(V\) is the windspeed is m/s and \(g\) is the gravitational constant.
- \( \hat{\mathbf{\omega}}\) is the normalised wind direction in the \(xz\) plane.

In simple terms, \( P(\mathbf{k}) \) can be thought of as the relative amplitude of the wave with wave vector \( \mathbf{k} \), based on the current windspeed and direction. The wave vector \( \mathbf{k} \) can take values in the range 

$$
\begin{gather}
\mathbf{k} = (k_x,\: k_z) = (2n\pi/L_x,\: 2m\pi/L_z) \\
n \in \mathbb{Z}, -N/2 \leq n < N/2 \\
m \in \mathbb{Z}, -M/2 \leq m < M/2 
\end{gather}
$$

In other words, \( \mathbf{k}\) is always a multiple of the fundamental frequency with respect to each axis. This results in our generated heightfield tiling perfectly, which can either be a drawback or a positive depending on the use case.

Going forward, we represent the frequency domain as an \(N \times M\) texture where each pixel corresponds to a wave vector, treating the centre of the image as the \(\mathbf{k} = \mathbf{0}\) vector. In terms of the integer pixel ids \(0 \leq p_x < N\) and \(0 \leq p_z < M\), we can calculate \(\mathbf{k}\) as 

$$
\mathbf{k} = (k_x, \: k_z) = (\frac{2\pi}{L_x}(p_x - \frac{N}{2}),\: \frac{2\pi}{L_z}(p_z - \frac{M}{2}))
$$

As an example, below is a Phillips Spectrum generated with the parameters \(C = 1\), \(V=20\), \(\mathbf{\omega} = (1, \:0.7)\), \(L_x = L_z = 300\) and \(N = M = 128\).
{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/phillips_128.png" alt="drawing" width="500"/>
            <figcaption> \( P(\mathbf{k}) \) represented in grayscale</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}
We can see that the spectrum is aligned with our chosen wind direction due to the dot product factor in \(P(\mathbf{k})\) eliminating the contribution of wave vectors perpendicular to \(\hat{\mathbf{\omega}} \). This effect can be emphasised by increasing the exponent of this factor, resulting in a spectrum (and therefore heightfield) more visibly aligned with \( \hat{\mathbf{\omega}} \).
{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/phillips_128_8.png" alt="drawing" width="500"/>
            <figcaption> \( P(\mathbf{k}) \) using a cosine factor of \( | \hat{\mathbf{k}} \cdot \hat{\mathbf{\omega}} |^8\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}


## Building the Initial Spectrum 
Armed with our wave spectrum, it would now be good to introduce some unpredictability to each wave, both in its phase and amplitude. We will achieve this by generating a random complex number \(\mathbf{z} = z_r + iz_i\) for each wave (or, equivalently, each pixel in our spectrum texture), with \(z_r\) and \(z_i\) being independent draws from a standard normal distribution with \(\mu = 0\) and \(\sigma = 1\). Tessendorf states that values generated in this way are a close fit for real experimental data on ocean waves. These values are generated only once on the startup of the program to act as the "seed" for the system. In terms of implementation, these values can be generated by feeding \(MN\) pairs of independent Uniform \(U[0,1]\) samples to a compute shader and apply the [Box-Muller](https://en.wikipedia.org/wiki/Box–Muller_transform) method to transform them into independent \(\mathcal{N}(0,1)\) samples.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/random_128.png" alt="drawing" width="500"/>
            <figcaption>A texture filled with random values, generated as detailed above. Here, the red channel represents the real component and the green channel represents the imaginary component.</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

This random data is then combined with the previously generated spectrum and a normalising factor to form our initial wave spectrum, denoted \(\tilde{h}_0(\mathbf{k})\) and calculated as follows:

$$
\tilde{h}_0(\mathbf{k}) = \frac{1}{\sqrt{2}}(z_r + iz_i)\sqrt{P(\mathbf{k})}
$$

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/InitialSpec_128.png" alt="drawing" width="500"/>
            <figcaption>A texture filled with random values, generated as detailed above. Here, the red channel represents the real component and the green channel represents the imaginary component.</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}


## Animating the Spectrum in Time
The next hurdle to overcome is how we take our static spectrum and progress it through time. While working in the frequency domain, moving forwards in time corresponds to changing the phase of each wave by a constant rate. As the phase of each wave is represented by its complex argument, we can achieve this by rotating each wave's complex amplitude at an appropriate speed (i.e. multiplying by some \(\mathbf{q} \in \mathbb{C}\) with \(| \mathbf{q} | = 1\)). To determine the speed at which we will progress each wave, we use a **dispersion relation** - a function that relates the wavelength, \(k\), of our wave-vector to its frequency. Following Tessendorf's implementation, we will use the dispersion relation 

$$
\omega(k) = \sqrt{gk}
$$

which is appropriate when dealing with ocean waves in deep water, as is the case here.

Using this dispersion relation, one way to obtain an animated spectrum \(\tilde{h}(\mathbf{k}, t)\) would be to take 

$$
\tilde{h}(\mathbf{k}, t) = \tilde{h}_0(\mathbf{k}) \exp (i\omega (k)t)
$$

Taking this as our animated spectrum would work fine, however there is one key issue with this approach. The issue is that this version of \(\tilde{h}(\mathbf{k}, t)\) is not hermitian for constant \(t\) i.e. we do not have that, in general, \(\tilde{h}(\mathbf{-k}, t) = \tilde{h}^*(\mathbf{k}, t)\). The result of this is that the (inverse )DFT of \(\tilde{h}\), \( \mathcal{F}^{-1}[\tilde{h}]\), is not necessarily a real-valued function. While we could just discard the non-real component after taking the DFT, it makes more intuitive sense for the heightfield to be a purely real-valued function. More importantly however, ensuring that the output of our FFTs are real-valued functions will pave the way for some powerful optimisations further down the line. To make the spectrum hermitian, and by extension the heightfield purely real, we modify it by also propagating waves in the opposite direction to obtain the final version of the spectrum:

$$
\tilde{h}(\mathbf{k}, t) = \tilde{h}_0(\mathbf{k}) \exp (i\omega (k)t) + \tilde{h}_0^*(-\mathbf{k}) \exp (-i\omega (k)t)
$$

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/movingSpec.gif" alt="drawing" width="500"/>
            <figcaption>\(\tilde{h}(\mathbf{k}, t)\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

For implementation purposes, it is convenient to pack both the initial spectrum \(\tilde{h}_0(\mathbf{k})\) and the conjugated initial spectrum \(\tilde{h}_0^*(\mathbf{-k})\) into one texture using all four channels. From this texture, the spectrum at any time \(t\) can be directly calculated using the expression above.
{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/initSpecPair.png" alt="drawing" width="500"/>
            <figcaption> \(\tilde{h}_0(\mathbf{k})\) and \(\tilde{h}_0^*(\mathbf{-k})\) packed using all four texture channels (note alpha is not visible)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}


## Moving into the Spatial Domain
With an animated hermitian spectrum, we are now ready to "pull out" our waves into the time domain and generate the heightfield texture. We will denote the height at a point \(\mathbf{x}\) and at time \(\mathbf{t}\) by \(h(\mathbf{x}, t)\), where \(\mathbf{x}\) takes on the discrete values 

$$
\begin{gather}
\mathbf{x} = (\frac{nL_x}{N}, \frac{mL_z}{M}) \\
n \in \mathbb{Z}, -N/2 \leq n < N/2 \\
m \in \mathbb{Z}, -M/2 \leq m < M/2
\end{gather}
$$

In terms of the pixel coordinates \((p_x, \: p_z)\) we can also express \(\mathbf{x}\) as 

$$
\mathbf{x} = (\frac{L_x}{N}(p_x - \frac{N}{2}), \frac{L_z}{M}(p_z - \frac{M}{2}))
$$

The function \(h(\mathbf{x}, t)\) is then calculated by taking a standard 2D DFT of the moving spectrum as follows:

$$
h(\mathbf{x}, t) = \sum_{\mathbf{k}} \tilde{h}(\mathbf{k}, t) \exp(i\mathbf{k} \cdot \mathbf{x})
$$

where we are summing over all wave vectors in the spectrum texture. For now, we can simply implement the above expression directly with a regular DFT (not FFT) using a compute shader. This is fast enough for smaller texture sizes, however if we wish to use texture resolutions of \(512\) or above for more realistic results, we must upgrade to the FFT.

{{< rawhtml >}}
    <center>
        <figure>
            <img class="pixelated" src="/images/fftWater_1/heightmap.gif" alt="drawing" width="500"/>
            <figcaption>\(h(\mathbf{x}, t)\) represented by the red channel with \(N=M=128\)</figcaption>
        </figure>
    </center>
{{< /rawhtml >}}

## References
- [Jerry Tessendorf: Simulating Ocean Water](https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf)
- [Thomas Gamper: Ocean Surface Generation and Rendering](https://www.cg.tuwien.ac.at/research/publications/2018/GAMPER-2018-OSG/GAMPER-2018-OSG-thesis.pdf)
- [Fabio Suriano: An introduction to Realistic Ocean Rendering through FFT (Talk)](https://www.youtube.com/watch?v=ClW3fo94KR4)
