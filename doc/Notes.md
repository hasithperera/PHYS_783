
## Basic notes from the code

- Energy data is exported using the function `energy_output()`
- Reading this data file from the example gives 36 columns
	- Consistant with the code in `vp_funcs.f90/energy_output` 

| variable name | usage                 | function/line no |
| ------------- | --------------------- | ---------------- |
| nspec         |                       |                  |
| etot          | total energy          |                  |
| ephi          | field energy          |                  |
| efs(1:nspec)  | total dist fun energy |                  |
| efsw          | wave dist fun energy  |                  |
|               |                       |                  |

---

## 04/12/2024

- Why is the code separated into phase space densitied
	- 
- How can I export phase phase densities

- Initialize data
	- `vp_data.f90` line 270 `read_in_params`
- output file formatting
	- phase_space.x :  contains distribution function data
		- $nx*nv$ : data points for one time step
---
- $\phi(x,t)$ data for stationary wave: run2 input file
- This is a standing wave. with time its not moving

![](res/Pasted%20image%2020240412122114.png)

---
- $n(x)$ for a propagating wave: run3 input file
![](res/Pasted%20image%2020240412155543.png)

---
## initialization of different species (3 and 4)

- Applying a bulk flow works
- Temperature changes are not working. (changing that doesn't change anything on the initialization end)
	-  The implementation looks like change in temperature in not implemented which is surprising 
	![](res/Pasted%20image%2020240422125138.png)
	- #ask_paul
- find $v_\phi$ in the vp code using the initial conditions and the notmalization

---
# initial results 

- `landaudamp4_0.in`: shows landau damping and how the energy changes in the system
	![500](res/Pasted%20image%2020240423114812.png)
---

- compare data from `landaudamp4.in` run. This run has a beam on the tail one one side of the wave
- comparing with `landaudamp4_1.in`
	![500](res/Pasted%20image%2020240423122737.png)
	- orange data if for the initial landau damping only run
	- Clearly the wave now take much longer to damp out (comparing these two rates would be one way to get an estimate of the growth rate)
	- Still this simulation has two waves left and right and only one beam that gives it a growth 
---
- `landaudamp4_2.in`: implement a total of 6 species to apply a bump on the $-v_\phi$ side too
	![350](res/Pasted%20image%2020240423125256.png)
	- Damping rate and the growth rate is almost balanced. Check the theory value for this and see it this is consistent

---
## chat with Paul

Explicit vs implicit numerical algorithms
- consider a diffusion equation for example $\frac{\partial f}{\partial t} = c \frac{\partial f}{\partial x}$
- We can write two finite difference schemes 
- explicit: If you use know steps
 $$\frac{f^{n+1}_j - f^{n}_j}{\Delta t}=c\frac{f^{n}_{j+1}-f^{n}_{j}}{2\Delta x}$$
 - implicit
 $$\frac{f^{n+1}_j - f^{n}_j}{\Delta t}=c\frac{f^{n+1}_{j+1}-f^{n+1}_{j}}{2\Delta x} $$

Implicit code will not crash even if you have a big time step. depending on what your need you can choose the correct time step. if you want to study reconnection often people use implicit codes with big time steps. specially when you want to study reconnection and after its come to a steady state.

---
# Propagating wave runs

- `landaudamp8`: using the initial conditions `ic=2`
	![350](res/Pasted%20image%2020240424093151.png)
	- looking at $\int f(x,v,t) dx$
	- clearly can see particles being sent to higher energy >> Landau damping signature
---
- `landaudamp8_1`: add an additional distribution
	![350](res/Pasted%20image%2020240424093645.png)
	- my species_1 is acting as the same. I think this is expected. species are evolved independently so it will not work as expected
---
## try to init a bump on s=0

- When I initialize the bump for s=0 it now shows the correct f as time goes on 
- still $\phi$ is not working (it still shows a decay)
	![width:550](res/Pasted%20image%2020240424130504.png)

---
## Evaluate the growth rate

$$
 \omega_i = \frac{\pi}{2} \omega_{pe} \frac{4\pi e^2}{k^2m_e} \frac{d F_{e0}}{d v_z}
$$

For the bump population
$$
f_{bump} = n e^{\left( -\frac{v^2}{2 v_{th}^2}\right)}
$$
- $v_{th}$ definition in the vp code uses a 2 outside

## Simulation details

### Basic landau damping runs

| Sim id | nt   | $v_{max}$ | $nv$ | $nx$ | $nl$ (sets kl_d) | $\alpha$ | vp                       | notes                                                  | ic  |
| ------ | ---- | --------- | ---- | ---- | ---------------- | -------- | ------------------------ | ------------------------------------------------------ | --- |
| 2      | 5000 |           |      |      |                  |          | org                      | initial week damping simulation with nonlinearity      | 0   |
| 2_1    |      |           |      |      |                  |          | mod by adding bump on +v | modify v binning to +/- 12 v_th<br>add f bump at v=7.5 |     |
| 3      | 5000 | 5         | 256  | 128  | 2                |          | org                      | moderate damping                                       | 1   |

---
### Testing single wave

| Sim id | nt    | $v_{max}$ | $nv$ | $nx$ | $nl$ (sets kl_d) | $\alpha$ | vp  | notes                                   | ic  |
| ------ | ----- | --------- | ---- | ---- | ---------------- | -------- | --- | --------------------------------------- | --- |
| 4      | 10000 | 10        | 256  | 128  | 2                | 4        |     |                                         | 1   |
| 4_0    |       |           |      |      |                  | 2        |     | amp of bump 0.5                         | 1   |
| 4_1    | 10000 |           |      |      |                  | 4        |     | amp of bump 1.0                         | 1   |
| 4_2    | 10000 |           |      |      | 2                | 6        |     | amp of bump = 0.1,4 extra species added | 1   |
| 4_3    |       |           |      |      |                  |          |     | same as 4_2. low amp bump               |     |

---
#### 4_0 results

![bg left width:500](res/Pasted%20image%2020240430094342.png)
- Stack plot of $\phi(x)$
- Decay in color indicates the wave perturbation dropping as time goes on 

---

- $\delta f$ integrated over space shows particle energization near the resonant frequency.
- lighter colors are early in time
- black dashed line is the resonant frequency

![bg left width:500](res/Pasted%20image%2020240430084046.png)

---
#### 4_1 results

- Energy plot shows slower decay. but $f$ looks odd. Moving on to implementing a bump on the $-3.5 v_{th}$ side to see if energy will be given back to the fields

---
#### 4_2 results

- Two bumps are implemented using 4 additional species
- ion and co-located electron bump to maintain charge neutrality 

![bg right width:600px](res/Pasted%20image%2020240430094754.png)

---
#### 4_3 results

- Same simulation as 4_3 but the bump density was dropped to $0.5$
- Energy exchange plot shows slow decay compared to ref Landau damping sim
  
![bg right width:600](res/Pasted%20image%2020240430084856.png)

---
## To-do

- [x] simulate a propagating wave. #task ✅ 2024-04-24
- [x] add a new population to the existing code #task ✅ 2024-04-24
- [ ] investigate if putting a new population vs same population impact #task
