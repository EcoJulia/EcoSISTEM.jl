# Model development

*Description of the model as it currently stands.*

## Epidemiology

* We currently have an SEI2HRD compartmental model: Susceptible-Exposed-Asymptomatic-Symptomatic-Hospitalised-Recovered-Dead, as seen in Figure 1 below.

### Virus update loop
* Virus grows as the sum of several Binomal draws from each infection category,
$NV_{k, c}$ ~ $Binomial(V_k, p_c)$
where $NV$ is new virus, $k$ is the current grid square, $c$ is the disease class, and $p_c$ is the probability of virus generation per disease class.
* The probability of generating virus is dependent on the population size of that disease class ($P_{c,k}$), the growth parameter of that disease class ($g_c$), and the match of the virus to the environment at that location ($T_k$):
$p_c = g_c * P_{c,k} * T_k$
* The virus decays similarly, according to a set probability $d$ and the inverse match with environment ($T_k^{-1}$):
Decayed virus ~ $Binomial(V_k, d * T_k^{-1})$
* The newly generated virus is distributed in space (i.e. grid square $k \rightarrow j$) via a Gaussian kernel:
$D_c(k,j)*NV_{k,c}$ where $D_c$ is a gaussian kernel per infection class.
* Currently, only the symptomatic and asymptomatic infectious classes can generate and spread virus, and $p_c = 0$ for all others.

### Disease class update loop
* Birth/death per class: Susceptibles are born into the population at a set probability per class. There is also a background probability of death from each disease class.
  * $Births_{k,c}$ ~ $Binomial(P_{k,c}, b_c)$
  * $Deaths_{k,c}$ ~ $Binomial(P_{k,c}, d_c)$
* Transitions per class: transitions between disease classes happen according to a transition matrix $M$, which are calculated as moves into the disease class from other categories:
$in_{k, c} = \sum_{j}^{J} Binomial(P_{k, j}, M_{c, j})$
 where $j \in{1,..,J}$ number of classes.
* $M$ is constructed in advance and is altered for movement from Susceptible to Exposed categories by the amount of virus available in that location, $V_k$.



![](SEI2HRD.png)
*Figure 1: Current model structure.*

## Code Structure

* The virus update loop happens first, parallelised over disease category. The virus must move between locations, so each process must have full access to the entire space. This will parallelise better when we have more disease categories, like age.
* The disease class update loop happens second, parallelised over space. There is no movement between locations, but instead transitions between different categories.
* An overall matrix housing the abundances per grid square and per disease category is housed in the `EpiSystem`, along with information on environment and a lookup table of moves between different grid squares for each kernel.
* At every iteration of the update step, this matrix is updated in place to avoid additional memory allocation.
