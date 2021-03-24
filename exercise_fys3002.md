---
title:
-	"Exercise in FYS-3002"
author:
-	Eirik Rolland Enger
tags:	[UiT]
date:
-	\today
geometry: margin=2.5cm
autoEqnLabels: true
numbersections: true
header-includes: |
	\usepackage{siunitx}
	\usepackage{cleveref}
output:
	pdf:document:
		templete: NULL

---

# Solving the equation for the ISR spectrum
> Following @Farley1999. Equations from @Farley1999 are written FH (xx).

We are going to derive the equation for the power spectral density of the electron number density,
namely the expression

$$
\langle |n_\mathrm{e}(\boldsymbol{k}, \omega)|^2\rangle=\frac{n_{\mathrm{e},0}}{\pi\omega}\frac{
\Im\{F_\mathrm{e}\}|1+\chi_\mathrm{i}|^2+\Im\{F_\mathrm{i}\}|\chi_\mathrm{e}|^2
}{
|1+\chi_\mathrm{e}+\chi_\mathrm{i}|^2
}.
$$ {#eq:electron_number_density}

## Vlasov's equation

Let us start from the Vlasov's equation.
This is similar to the Boltzmann equation, but we omit the collision term.
The Vlasov equation can be written as (FH 4.24)
$$
\partial_t f_\alpha + \boldsymbol{v}\cdot \partial_{\boldsymbol{r}}f_\alpha+\mu_\alpha
[\boldsymbol{E}+\boldsymbol{v}\times\boldsymbol{B}]\cdot\partial_{\boldsymbol{v}}f_\alpha=0
$$
where $f=f(\boldsymbol{r},\boldsymbol{v},t)$ describe the phase space, $\mu_\alpha$ is the
charge-to-mass ratio of particle species $\alpha$ and $\boldsymbol{E}$ and $\boldsymbol{B}$ are the
electric and magnetic fields, both functions of space and time.

Let us assume all parameters to consist of a linear term and a higher order term, that is, we assume
all parameters are on the form $f=f_0[1+f_1]$ where $f_0$ is linear and $f_1$ is non-linear and that
$f_1 \ll 1$. Let us also write up the Fourier transform and Laplace transform of $f_1$:
$$ \begin{aligned} f_1(\boldsymbol{r},\boldsymbol{v},t)&=\sum_{\boldsymbol{k}}
f_1(\boldsymbol{k},\boldsymbol{v},t)\exp(-i\boldsymbol{k}\cdot\boldsymbol{r})\\
f_1(\boldsymbol{k},\boldsymbol{v},s)&=\int_0^\infty
f_1(\boldsymbol{k},\boldsymbol{v},t)\exp(-st)\mathrm{d}t. \end{aligned} $$ We linearize the Vlasov
equation and obtain $$ sf_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s)
-f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},t=0) -i\boldsymbol{k} \cdot
\boldsymbol{v}f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s) +\mu_\alpha \left[
\frac{1}{f_{\alpha,0}(\boldsymbol{v})}\boldsymbol{E}
\cdot\partial_{\boldsymbol{v}}f_{\alpha,0}(\boldsymbol{v}) -\boldsymbol{B}
[\boldsymbol{v}\times\partial_{\boldsymbol{v}}f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s)]
\right]=0. $$ It can be shown that this has solution $$
f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s) =\frac{1}{\mu_\alpha B} \int_{-\infty}^\varphi
g_\alpha (\varphi, \varphi') \left\{ f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v}',t=0)\mp
\frac{i2X_\mathrm{p}^2}{f_{\alpha,0}(\boldsymbol{v}')}\boldsymbol{k} \cdot\boldsymbol{v}'
[Zn_\mathrm{i}(\boldsymbol{k},s) - n_\mathrm{e}(\boldsymbol{k},s)] \right\} \mathrm{d}\varphi'
$${#eq:perturbation} where $Z$ is the ion charge in units of the elementary charge. The primes (e.g.
on $\boldsymbol{v}'$) refer to terms on the *unperturbed* orbit. Specifically we have that the
unperturbed velocity is $$ \boldsymbol{v}' = \mathbf{e}_1 w\cos\varphi'+\mathbf{e}_2 w\sin\varphi'
+\mathbf{e}_3 u $$ and we see that the velocity is a function of the variable $\varphi$ (i.e.,
$f_1(\boldsymbol{k},\boldsymbol{v},s)=f_1(\boldsymbol{k},w,u,\varphi,s)$). $g(\varphi,\varphi')$ can
be seen as an integrating factor [@Bernstein1958].

Density perturbations can then be obtained by integration:
$$
n_\alpha(\boldsymbol{k},s)=\int
f_{\alpha,0}(\boldsymbol{v}')f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s)\mathrm{d}\boldsymbol{v}
=n_{\alpha,0}\int f_{\alpha,1}(\boldsymbol{k},\boldsymbol{v},s)\mathrm{d}\boldsymbol{v}.
$$

## Electron number density

The function $F_\alpha$ is defined as
$$
\begin{aligned} F_\alpha(\boldsymbol{k},\omega) &=1+i\omega\int_0^\infty \exp\left\{ i\omega
\tau-\frac{k_\mathrm{r}^2 T_\alpha k_\mathrm{B}\sin^2\theta}{m_\alpha\Omega_\alpha^2}[1-\cos(\omega
\tau)]-\frac{1}{2}(k_\mathrm{r}\tau\cos\theta)^2\frac{T_\alpha k_\mathrm{B}}{m_\alpha}
\right\}\mathrm{d}\tau\\ &=1+i\omega G_\alpha(\boldsymbol{k},\omega). \end{aligned}
$$
The integral $G_\alpha(\boldsymbol{k},\omega)$ is often referred to as a Gordeyev integral. Another
useful parameter is: $$
\chi_\alpha(\boldsymbol{k},\omega)=\frac{1}{k^2\lambda_\mathrm{D}^2}F_\alpha(\boldsymbol{k},\omega)
=\frac{1}{k^2\lambda_\mathrm{D}^2}[1+i\omega G_\alpha(\boldsymbol{k},\omega)]. $${#eq:chi}

Using our definitions of $F_\alpha$ and $\chi_\alpha$ we can now readily compute the IS spectrum in
eq. (@eq:electron_number_density): $$ \langle |n_\mathrm{e}(\boldsymbol{k}, \omega)|^2\rangle
=\frac{n_{\mathrm{e},0}}{\pi\omega}\frac{
\Im\{F_\mathrm{e}\}|1+\chi_\mathrm{i}|^2+\Im\{F_\mathrm{i}\}|\chi_\mathrm{e}|^2 }{
|1+\chi_\mathrm{e}+\chi_\mathrm{i}|^2 }. $${#eq:number_density_2} We may also rewrite the above
using that $\Im\{F_\alpha\}=\omega G_\alpha$: $$ \langle |n_\mathrm{e}(\boldsymbol{k},
\omega)|^2\rangle =\frac{n_{\mathrm{e},0}}{\pi}\frac{
G_\mathrm{e}|1+\chi_\mathrm{i}|^2+G_\mathrm{i}|\chi_\mathrm{e}|^2 }{
|1+\chi_\mathrm{e}+\chi_\mathrm{i}|^2 }=\frac{n_{\mathrm{e},0}}{\pi}\left[\frac{
G_\mathrm{e}|1+\chi_\mathrm{i}|^2 }{ |1+\chi_\mathrm{e}+\chi_\mathrm{i}|^2 }+\frac{
G_\mathrm{i}|\chi_\mathrm{e}|^2 }{ |1+\chi_\mathrm{e}+\chi_\mathrm{i}|^2 } \right].
$${#eq:number_density_3}

# Problems
## *Problem 1* {-}
### Task 1 {-}
What is the physical interpretation of $\chi$ in eq.\ (@eq:chi)?

### Task 2 {-}
Explain what the term $f_{\alpha, 1}$ in eq.\ (@eq:perturbation) describe.

### Task 3 {-}
Explain what the two terms in the square bracket in eq.\ (@eq:number_density_3) describe.

## *Problem 2* {#prob:2 -}

Using any one of the expressions for $\langle |n_\mathrm{e}(\boldsymbol{k}, \omega)|^2\rangle$,
write a program that calculates the power spectral density. The program should accept a number of
input parameters:

> $f_\mathrm{r}$
> : Radar frequency
>
> $n_\mathrm{e}$
> : Electron number density
>
> $B$
> : Magnetic field strength
>
> $m_\mathrm{i}$
> : Ion mass
>
> $T_\mathrm{e}$
> : Electron temperature
>
> $T_\mathrm{i}$
> : Ion temperature
>
> $\theta$
> : Aspect angle (the angle between the radar pointing direction and the magnetic field)

Explain the code. Some things to consider:

-   Where in the code were the different equations solved?
-   How was the numerical calculation implemented? 

The code itself should be well commented and included as an appendix.

*Hint*: Before you integrate all the way to infinity: where along the axis of integration does most
of the information lie?

## *Problem 3* {#prob:3 -}

We will now look at some specific parameters using our program. Run your program with the parameters
given as:

 Parameter     |  Unit            |  Value
:---------:    | :----:           | :-----:
$f_\mathrm{r}$ | [$\si{\hertz}$]  | $\num{430e6}$
$n_\mathrm{e}$ | [$\si{m^{-3}}$]  | $\num{2e10}$
$B$            | [$\si{\tesla}$]  | $\num{3.5e-5}$
$m_\mathrm{i}$ | [$\mathrm{amu}$] | $16$
$T_\mathrm{e}$ | [$\si{\kelvin}$] | $200$
$T_\mathrm{i}$ | [$\si{\kelvin}$] | $200$
$\theta$       | [$\si{\degree}$] | $135$
for frequencies $f\in[-\num{2e6},\num{2e6}]$.

### Task 1 {-}

Where could an experiment with these parameters be done? Make a sketch that include the radar beam
and the magnetic field line of Earth, in addition to the radar's approximate position on earth.
Assume that the radar points directly upwards, i.e., towards zenith.

### Task 2 {-}

The spectrum is plotted for frequencies $f\in[-\num{2e6},\num{2e6}]$; relative to an observer at the
radar location, which way does the features found at positive frequencies in the spectrum move?

### Task 3 {-}

Plot the resulting power spectra calculated by the program and explain what the different peaks
represent.

## *Problem 4* {#prob4 -}
In this exercise we will use the parameters:

 Parameter     |  Unit            |  Value
:---------:    | :----:           | :-----:
$f_\mathrm{r}$ | [$\si{\hertz}$]  | $\num{933e6}$
$n_\mathrm{e}$ | [$\si{m^{-3}}$]  | $\num{2e11}$
$B$            | [$\si{\tesla}$]  | $\num{5e-5}$
$m_\mathrm{i}$ | [$\mathrm{amu}$] | $16$
$T_\mathrm{i}$ | [$\si{\kelvin}$] | $2000$
$\theta$       | [$\si{\degree}$] | $180$

### Task 1 {-}

Calculate the power spectral density on $f\in[\num{3.5e6},\num{7e6}]$ for
$T_\mathrm{e}=\SI{2000}{\kelvin},\,\SI{4000}{\kelvin},\,\SI{6000}{\kelvin},$ and
$\SI{8000}{\kelvin}$ and plot the power spectra.

### Task 2 {-} 

Explain the changes that can be seen as the electron temperature change.

### Task 3 {-}

Explain what the assumption $k^2\lambda_\mathrm{D}^2\ll 1$ is referring to. Is this assumption valid
for all temperatures? Why/why not?

## *Problem 5* {#prob5 -}
You should now be able to experiment a bit for yourself.

### Task 1 {-}

Explain which parameter(s) that needs to be changed to obtain a similar plot as is shown in
\cref{fig:ionlines}.

### Task 2 {-}

Reproduce the plot in \cref{fig:ionlines} using your own program. Include values on the axis and
labels to all spectra.

![Power spectral density plot \label{fig:ionlines}](ionline.pdf){ width=80% }

## References

