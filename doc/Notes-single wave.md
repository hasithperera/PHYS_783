### Single traveling wave : Simulation details (internal use)

| run_id | nt    | cfl   |     |     |     | $\alpha$ |     | notes                       | ic  | eig_opt |
| ------ | ----- | ----- | --- | --- | --- | -------- | --- | --------------------------- | --- | ------- |
| 8      | 10000 |       |     |     |     |          |     | bump in s=1                 | 2   | 0       |
| 8_0    | 10000 | 0.125 |     |     |     | 2        |     | pure run                    | 2   | 0       |
| 8_1    | 10000 | 0.125 |     |     |     | 4        |     | single wave shows damping   | 2   | 0       |
| 8_2    | 9     |       |     |     |     | 4        |     | high amplitude bump on tail |     |         |
| 8_3    | 10000 | 0.05  |     |     |     | 2        |     | run from documentation      |     |         |
|        |       |       |     |     |     |          |     |                             |     |         |
|        |       |       |     |     |     |          |     |                             |     |         |
|        |       |       |     |     |     |          |     |                             |     |         |

---
## Landau damping only run 

- $n(x,t)$ looks good shows a traveling wave (correct periodicity)
- $\phi(x,t)$ should have a periodic boundary but that is not the case in the data
	- looked around at the green's function implementation. Don't fully understand the numerical algorithm there
![bg right width:450](res/Pasted%20image%2020240501110550.png)

---

## Energy plots

- Not great. This is a run from the documentation 
- $\tilde{n}=0.02$

![bg right width:650](res/Pasted%20image%2020240501110849.png)

---
## f(t)

- lighter colors early in time 
- right plot: $\int f(t=\infty,x) dx$

![](res/Pasted%20image%2020240501111157.png)

![bg right width:650](res/Pasted%20image%2020240501111603.png)