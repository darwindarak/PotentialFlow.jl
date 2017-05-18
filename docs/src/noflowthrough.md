# Enforcing No-Flow-Through

!!! warning
    Under construction...

```math
\def\ddt#1{\frac{\mathrm{d}#1}{\mathrm{d}t}}

\renewcommand{\vec}{\boldsymbol}
\newcommand{\uvec}[1]{\vec{\hat{#1}}}
\newcommand{\utangent}{\uvec{\tau}}
\newcommand{\unormal}{\uvec{n}}

\renewcommand{\d}{\,\mathrm{d}}

\newcommand{\cross}{\times}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\im}{\mathrm{i}}
\newcommand{\eu}{\mathrm{e}}
\newcommand{\pint}{\int}
\newcommand{\conj}[1]{#1^\star}
\newcommand{\Res}[2]{\mathrm{Res}\left(#1,#2\right)}
\newcommand{\real}[1]{\mathrm{Re}\left\{#1\right\}}
\newcommand{\imag}[1]{\mathrm{Im}\left\{#1\right\}}
```
We are interested in enforcing the no-flow-through condition on an infinitely thin, flat plate undergoing rigid body motion.
The plate can be parameterized by its length, ``L``, centroid position, ``\vec{c}``, and its angle of attack, ``\alpha``.
Its motion is then specified by its centroid velocity, ``\dot{\vec{c}}``, and angular velocity, ``\dot{\vec{\alpha}}``.

## Vortex Sheet Strength

The plate is represented with a bound vortex sheet that constantly adjusts its circulation to enforce no-flow-through on its surface.
We can show that the distribution of circulation, ``\gamma``, is governed by the following integral equation:
```@raw html
<details>
<summary></summary>
The no-flow-through condition requires that the component of fluid velocity normal to the sheet must be equal to the normal velocity of the sheet itself, i.e.
$$
\begin{align*}
    \unormal \cdot \vec{u}(\vec{x}_s)
& = \unormal \cdot \left[ \dot{\vec{c}} + \dot{\alpha} \cross (\vec{x}_s - \vec{c}) \right] \\
& = \left(\unormal \cdot \vec{c}\right) + \dot{\alpha} l
\end{align*}
$$
where

- $\vec{u}$ is the fluid velocity
- $\vec{x}_s$ is a position on the plate
- $\unormal$ is a unit vector normal to the plate
- $l \in [ -L/2, L/2 ] $ is distance between $\vec{x}_s$ from the plate centroid

We can decompose the velocity field at any point in the fluid into contributions from the bound vortex sheet, $\vec{u}_s$, and the free vorticity in the ambient fluid, $\vec{u}_A$:
$$
\vec{u}(\vec{x}) = \vec{u}_s(\vec{x}) + \vec{u}_A(\vec{x}),
$$
so the no-flow-through condition can be written as:
$$
\unormal \cdot \vec{u}_s(\vec{x}) = \left(\unormal \cdot \vec{c}\right) + \dot{\alpha} l - \unormal \cdot \vec{u}_A(\vec{x}).
$$

The velocity field induced by a vortex sheet, $\vec{u}_x(\vec{x})$, is given by
$$
\vec{u}_s(\vec{x}) = \frac{1}{2\pi}
\int_\mathcal{C} \gamma(l) \,\uvec{k} \cross
\frac{\vec{x} - \vec{x}_s(l)}{\abs{\vec{x} - \vec{x}_s(l)}^2}
\d{l}
$$
where

- $\gamma$ is the strength of the sheet
- $\mathcal{C}$ is the curve occupied by the sheet
- $\uvec{k}$ is the unit vector point out of the plane.

The position along the vortex sheet can be expressed as
$$
\vec{x}_s(l) = \vec{c} + l\utangent
$$
where $\utangent$ is the unit tangent along the sheet.
Similarly, since we are interested in evaluating the velocity along the sheet, we can write
$$
\vec{x}(l) = \vec{c} + \lambda\utangent.
$$
We can then write self-induced velocity of the bound vortex sheet as
$$
\vec{u}_s(\lambda) = \frac{\unormal}{2\pi}
\int_{-\frac{L}{2}}^\frac{L}{2} \frac{\gamma(l)}{\lambda - l}
\d{l}.
$$
Substituting this expression back into the no-flow-through condition, we get
</p>
</details>
```
```math
\begin{equation}
\frac{1}{2\pi}
\int_{-L/2}^{L/2} \frac{\gamma(\lambda)}{l - \lambda}
\d{\lambda}
= \unormal \cdot \vec{\dot{c}}
+ \dot{\alpha} l
- \unormal \cdot \vec{u}_A(l)
\label{eq:integral-equation}
\end{equation}
```
The solution to this integral equation can be found in [^Muskhelishvili].
If the velocity induced by ambient vorticity on the plate can be expanded into a Chebyshev series:
```math
\unormal \cdot \vec{u}_A[l(s)] = \sum_{n = 0} A_n T_n(s),
```
and ``\Gamma_A`` is the total circulation in the ambient fluid, then the solution to ``\eqref{eq:integral-equation}`` can be written as:
```@raw html
<details>
<summary> </summary>
To make it easier to work with Chebyshev series, we will apply a change of variables $s := \frac{2l}{L}$ so that the integral above goes from $-1$ to $1$:
$$
\frac{1}{2\pi}
\int_{-1}^1 \frac{\gamma(s)}{\sigma - s}
\d{s}
= \unormal \cdot \vec{\dot{c}}
+ \frac{\dot{\alpha}L}{2} \sigma
- \unormal \cdot \vec{u}_A(\sigma)
$$
From <a href="#footnote-Muskhelishvili">[Muskhelishvili]</a>, we have that if
$$
\frac{1}{\pi\im} \int \frac{\varphi(t)}{t - t_0} \d{t} = f(t_0)
$$
then
$$
\varphi(t_0) = \frac{1}{\pi\im\sqrt{t_0 - 1}\sqrt{t_0 + 1}}
\int \frac{\sqrt{t - 1}\sqrt{t + 1}}{t - t_0} f(t) \d{t}
+
\frac{P(t_0)}{\sqrt{t_0 - 1}\sqrt{t_0 + 1}}
$$
where $P$ is an arbitrary polynomial that must be chosen to satisfy far-field boundary conditions.

In our case, we have $\varphi := \im \gamma$ and
$$
f := 2\sum_{n = 0}^\infty A_n T_n(\sigma) - 2\unormal \cdot \vec{\dot{c}} - \dot{\alpha}L \sigma
$$
so
$$
\gamma(\sigma)
=
\frac{-2}{\pi\sqrt{1 - \sigma}\sqrt{1 + \sigma}}
\int_{-1}^1 \frac{\sqrt{1 - s}\sqrt{1 + s}}{s - \sigma}
\left(
\sum_{n = 0}^\infty A_n T_n(s) - \unormal \cdot \vec{\dot{c}} - \frac{\dot{\alpha}L}{2} s
\right) \d{s}
+
\frac{P(t_0)}{\sqrt{1 - \sigma}\sqrt{1 + \sigma}}
$$

The integral above is made of terms with the form
$$
\pint_{-1}^1
\frac{\sqrt{1 - s}\sqrt{1 + s}}{s - \sigma} T_n(s)
\d{s}
$$
which we can simplify using the properties of Chebyshev polynomials into
$$
\pint_{-1}^1
\frac{\sqrt{1 - s}\sqrt{1 + s}}{s - \sigma} T_n(s)
\d{s}
=
\begin{cases}
-\pi T_1(\sigma) & n = 0 \\
-\frac{\pi}{2} T_2(\sigma) & n = 1 \\
-\frac{\pi}{2} \left[T_{n+1}(\sigma) - T_{n-1}(\sigma)\right] & n \ge 2
\end{cases}.
$$
This gives us
$$
\gamma(\sigma)
=
\frac{-2}{\pi\sqrt{1 - \sigma}\sqrt{1 + \sigma}}
\left\{
-\pi A_0 \sigma
-\frac{\pi}{2} A_1
+\sum_{n = 1}^\infty -\frac{\pi}{2}A_n \left[T_{n+1}(\sigma) - T_{n-1}(\sigma)\right]
+ \pi \left(\unormal \cdot \vec{\dot{c}}\right)\sigma
+ \frac{\pi}{2}T_2(\sigma)\frac{\dot{\alpha}L}{2}
\right\}
+
\frac{P(t_0)}{\sqrt{1 - \sigma}\sqrt{1 + \sigma}}.
$$

We can find $P$ by satisfying Kelvin's circulation theorem.
This means that the amount of circulation contained in the bound vortex sheet should the negative of the circulation contained in the ambient vorticity, i.e.
$$
\Gamma_s := \int_{-\frac{L}{2}}^{\frac{L}{2}} \gamma \d{l} = -\Gamma_A
$$

Again, we use properties of Chebyshev polynomials to reduce the integral to
$$
\begin{align*}
\frac{L}{2}\int_{-1}^1 \frac{P(s)}{\sqrt{1 - s}\sqrt{1 + s}} \d{s} & = -\Gamma_A,
\end{align*}
$$
which means that
$$
P = -\frac{2\Gamma_A}{L\pi}.
$$

So the final expression for the bound circulation is:
</details>
```
```math
\begin{equation}
\gamma[l(s)] =
\frac{-\frac{2\Gamma_A}{L\pi} + 2(A_0 - \unormal \cdot \vec{\dot{c}}) T_1(s) + (A_1 - \frac{\dot{\alpha}L}{2})T_2(s)}{\sqrt{1 - s^2}} - 2\sqrt{1 - s^2}\sum_{n = 2}^\infty A_n U_{n-1}(s)
\label{eq:gamma}
\end{equation}
```
!!! note
    This might look more similar to results from thin-airfoil theory if we rewrite the Chebyshev polynomials using trigonometric functions:
    ```math
    \gamma[l(\theta)] =
    \frac{-\frac{2\Gamma_A}{L\pi} + 2(A_0 - \unormal \cdot \vec{\dot{c}}) \cos\theta + (A_1 - \frac{\dot{\alpha}L}{2})\cos(2\theta)}{\sin\theta} - 2\sum_{n = 2}^\infty A_n \sin(n\theta).
    ```
    The key difference is that we are free to relax the Kutta condition at the trailing edge.

## Circulation

In addition to the distribution of circulation along the plate, it will be useful to know the amount circulation contained between one end of the plate to an arbitrary point on its surface.
By definition, we have
```math
\begin{align*}
\Gamma(l) & = \int_{-L/2}^l \gamma(\lambda) \d{\lambda} \\
\Gamma[l(s)] &= \frac{L}{2}\int_{-1}^s \gamma[l(\sigma)] \d{\sigma}.
\end{align*}
```
We can integrate ``\gamma`` term by term to obtain:
```@raw html
<details>
<summary></summary>
In equation $\eqref{eq:gamma}$, the Chebyshev polynomial of the second kind in ther summation can be written in terms of Chebyshev polynomials of the first kind:
$$
2\sqrt{1 - s^2}U_{n-1}(s)  = \frac{T_{n-1}(s) - T_{n+1}(s)}{\sqrt{1 - s^2}}.
$$
This means that all the terms in equation $\eqref{eq:gamma}$ can be expressed in the form:
$$
\frac{T_n(s)}{\sqrt{1 - s^2}}.
$$
The integral of these terms are:
$$
\begin{align*}
\int_{-1}^s \frac{T_n(s)}{\sqrt{1 - s^2}} \d{s}
& = \int_{\cos^{-1} s}^\pi \cos(n\theta) \d{\theta} \\
& = \begin{cases}
\pi - \cos^{-1}s &: n = 0 \\
-\frac{1}{n}\sin\left(n\cos^{-1}s\right) &: n > 0
\end{cases}.
\end{align*}
$$
We can then multiply the expressions above with their corresponding coefficients to obtain:
</details>
```
```math
\Gamma[l(s)]
=\Gamma_A\left(\frac{\cos^{-1}s}{\pi} - 1\right) - \frac{L\sqrt{1 - s^2}}{2}\left[
2\left(A_0 - \unormal \cdot \vec{\dot{c}}\right)
+\left(A_1 - \frac{\dot{\alpha}L}{2}\right)s
+\sum_{n=2}^\infty
A_n \left(\frac{U_n(s)}{n+1} - \frac{U_{n-2}(s)}{n-1}\right)\right]
```

[^Muskhelishvili]: Muskhelishvili, Nikolaĭ Ivanovich, and Jens Rainer Maria Radok. Singular integral equations: boundary problems of function theory and their application to mathematical physics. Courier Corporation, 2008.
