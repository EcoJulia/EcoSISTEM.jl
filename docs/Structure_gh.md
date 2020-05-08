# Model development

*Description of the model as it currently stands.*

## Epidemiology

* We currently have an SEI2HRD compartmental model: Susceptible-Exposed-Asymptomatic-Symptomatic-Hospitalised-Recovered-Dead, as seen in Figure 1 below.

### Virus update loop
* Virus grows as the sum of several Binomal draws from each infection category, <img src="svgs/e377ad7e80565ddab0061ecf36275ff6.svg?invert_in_darkmode" align=middle width=41.486115pt height=22.381919999999983pt/> ~ <img src="svgs/7b3030c8cf6c610db383e71bf51e49de.svg?invert_in_darkmode" align=middle width=87.173295pt height=24.56552999999997pt/> where <img src="svgs/60491d9986c130c5330721a64b127048.svg?invert_in_darkmode" align=middle width=28.143225000000005pt height=22.381919999999983pt/> is new virus, <img src="svgs/63bb9849783d01d91403bc9a5fea12a2.svg?invert_in_darkmode" align=middle width=9.041505pt height=22.745910000000016pt/> is the current grid square, <img src="svgs/3e18a4a28fdee1744e5e3f79d13b9ff6.svg?invert_in_darkmode" align=middle width=7.087278000000003pt height=14.102549999999994pt/> is the disease class, and <img src="svgs/b17e856e76ef58f7655e6ace49d21d06.svg?invert_in_darkmode" align=middle width=14.092485000000003pt height=14.102549999999994pt/> is the probability of virus generation per disease class.
* The probability of generating virus is dependent on the population size of that disease class (<img src="svgs/7b28106ce905937c70c7bc86cdfc5386.svg?invert_in_darkmode" align=middle width=27.495435pt height=22.381919999999983pt/>), the growth parameter of that disease class (<img src="svgs/eb30b7762aec8ad6b9c53cf51888ac16.svg?invert_in_darkmode" align=middle width=13.663980000000002pt height=14.102549999999994pt/>), and the match of the virus to the environment at that location (<img src="svgs/345c3b641bdaa5a8313143aea23a0ecc.svg?invert_in_darkmode" align=middle width=16.809210000000004pt height=22.381919999999983pt/>): <img src="svgs/eec4544e03d5062f33136d895ddd6a89.svg?invert_in_darkmode" align=middle width=127.49269499999998pt height=22.381919999999983pt/>
* The virus decays similarly, according to a set probability <img src="svgs/2103f85b8b1477f430fc407cad462224.svg?invert_in_darkmode" align=middle width=8.524065000000002pt height=22.745910000000016pt/> and the inverse match with environment (<img src="svgs/ecaa4366e01bf34ec3aa9fbfc92b6749.svg?invert_in_darkmode" align=middle width=28.617105000000002pt height=28.839689999999997pt/>): Decayed virus ~ <img src="svgs/eea37afd8c2b8d0f226dc2ad4469ca0d.svg?invert_in_darkmode" align=middle width=161.733495pt height=28.839689999999997pt/>
* The newly generated virus is distributed in space (i.e. grid square <img src="svgs/b5ad87070466e0d57cbd063852c46855.svg?invert_in_darkmode" align=middle width=42.234225pt height=22.745910000000016pt/>) via a Gaussian kernel: <img src="svgs/1530e9afae284012a0f08ee808c4045a.svg?invert_in_darkmode" align=middle width=113.99058pt height=24.56552999999997pt/> where <img src="svgs/6e24f7e524feccaf8ee5c1016f0854bd.svg?invert_in_darkmode" align=middle width=19.411425pt height=22.381919999999983pt/> is a gaussian kernel per infection class.
* Currently, only the symptomatic and asymptomatic infectious classes can generate and spread virus, and <img src="svgs/200efead2cd5006b73e4cc6e345dbc79.svg?invert_in_darkmode" align=middle width=44.994675pt height=21.10812pt/> for all others.
* Newly generated virus infects other susceptibles either as an instantaneous force of infection (e.g. aerosol transmission) or through environmental transmission (e.g. through contact with a surface), both at separate rates. See below for more details.

### Disease class update loop
* Birth/death per class: Susceptibles are born into the population at a set probability per class. There is also a background probability of death from each disease class.
  * <img src="svgs/116b2fd65e6e5c80263945ae2a637458.svg?invert_in_darkmode" align=middle width=66.74266499999999pt height=22.745910000000016pt/> ~ <img src="svgs/92f895c4e75485aaa72ca8cc64fbd7f4.svg?invert_in_darkmode" align=middle width=132.67765500000002pt height=24.56552999999997pt/>
  * <img src="svgs/07e3266d0aef26caedb51aba291fe757.svg?invert_in_darkmode" align=middle width=70.30617pt height=22.745910000000016pt/> ~ <img src="svgs/418b3d7001f48d6981368ac51dfe3979.svg?invert_in_darkmode" align=middle width=134.173215pt height=24.56552999999997pt/>
* Transitions per class: transitions between disease classes happen according to a transition matrix <img src="svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.67348pt height=22.381919999999983pt/>, which are calculated as moves into the disease class from other categories: <img src="svgs/5bb225c2a9a7cd676904f753d364d995.svg?invert_in_darkmode" align=middle width=236.25574499999996pt height=32.19743999999999pt/> where <img src="svgs/f74e7a22d321a4aea27f9f75bc505885.svg?invert_in_darkmode" align=middle width=70.256175pt height=22.381919999999983pt/> number of classes.
* <img src="svgs/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode" align=middle width=17.67348pt height=22.381919999999983pt/> is constructed in advance and is altered for movement from Susceptible to Exposed categories by the amount of virus available in that location, <img src="svgs/f3d0ba8ca71eac08d173274ee914a043.svg?invert_in_darkmode" align=middle width=16.792215000000002pt height=22.381919999999983pt/>, or the instantaneous force of infection, <img src="svgs/a178eaacca2a06509317f53ac22ba70a.svg?invert_in_darkmode" align=middle width=17.770170000000004pt height=22.381919999999983pt/>, which disappears at the end of each step.



![](SEI2HRD.svg)
*Figure 1: Current model structure.*

## Code Structure

* The virus update loop happens first, parallelised over disease category. The virus must move between locations, so each process must have full access to the entire space. This will parallelise better when we have more disease categories, like age.
* The disease class update loop happens second, parallelised over space. There is no movement between locations, but instead transitions between different categories.
* An overall matrix housing the abundances per grid square and per disease category is housed in the `EpiSystem`, along with information on environment and a lookup table of moves between different grid squares for each kernel.
* At every iteration of the update step, this matrix is updated in place to avoid additional memory allocation.
