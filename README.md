# How to calculate the probability of lethality of a whale-ship collision using *whalestrike*
This tutorial aims to demonstrate how to use the R package [whalestrike](https://joss.theoj.org/papers/10.21105/joss.06473) for risk assessment of vessel strikes, particularly using AIS data.

## whalestrike R package
The R package [whalestrike](https://joss.theoj.org/papers/10.21105/joss.06473) was created to support the model presented by [Kelley et al. (2021)](https://onlinelibrary.wiley.com/doi/abs/10.1111/mms.12745). The model simulates a collision between an adult North Atlantic right whale (*Eubalaena glacialis*) and a ship and calculates the maximum reactive force resulting from that collision, which is then converted to an estimate of the probability of lethal injuries. Parameters such as mass, impact area, and speed of the ships are required as input for the ship dimensions and speed at the time of the collision. A Shiny app (`app2()`) can be used to view the effect of the different parameters on the lethality.

### Installation
Details on the package installation and functions are explained in [Kelley 2024](https://joss.theoj.org/papers/10.21105/joss.06473) and in the package [vignette](https://dankelley.github.io/whalestrike/).
```r
library(remotes)
install_github("dankelley/whalestrike", ref="main")

library(whalestrike)
```

### Functions
The probability of lethality is expressed as a **lethality index** over the time of the collision. The index is measured from the stress of the whale compression.  
  
`lethalityIndexFromStress(stress)`  
 
where `stress` is the whale compression stress in Pascals. The values of this logistic model can be found in the results of a strike simulation (`strike` function).  
  
The `strike` function simulates a collision between a vessel with specific characteristics and a whale with specific characteristics. The model inputs include a vector of time, the state and parameters of the vessel and the whale.  
  
`strike(t, state, parms)`  
   
| Arguments | Description |
|-|--------|
|t|vector of time (s)|  
|state|list containing:| 
|*xs*|*ship position (m)*|
|*vs*|*ship speed (m/s)*|
|*xw*|*whale position (m)*|
|*vw*|*whale speed (m/s)*| 
|parms|list containing the model parameters (created with `parameters`)|  
  
  
The function `parameters` controls the parameters of the vessel and whale characteristics in the simulation.  
  
  
`parameters(ms, Ss, Ly, Lz, species, lw, mw, Sw, l, a , b, s, Cs, Cw, logistic)`
  
| Arguments | Value | Description |
|--|----|--------|
|ms|*variable*, default = 45000|ship mass (kg)|
|Ss|*variable*|ship wetted area (m<sup>2</sup>) (if not given, it will be estimated from the ship mass)|
|Ly|*variable*, default = 1.15|ship impact horizontal extent (m)|
|Lz|*variable*, default = 1.15|ship impact vertical extent (m)|
|species|"N. Atl. Right Whale"|whale species|
|lw|13.7 m|whale length (m)|
|mw|NULL|whale mass (kg) (if not given, it will be estimated from the whale length)|
|Sw|NULL|whale surface area (m<sup>2</sup>) (if not given, it will be calculated from the whale length|
|l|default = c(0.025, 0.16, 1.12, 0.1)|vector of 4 of the thickness (m) of skin, blubber, sublayer and bone|
|a, b|default = 1e6 * c(19.6, 0.255, 0.255, 22.9)|vector of 4 giving values of the stress-strain law ($stress = a*(exp(b*strain)-1)$)|
|s|default = 1e6 * c(19.6, 0.255, 0.255, 22.9)|vector of 4 giving values of ultimate strength (Pa)|
|theta|default = 55|whale skin deformation angle (deg)|
|Cs|default = 1e-2|drag coefficient for ship (calculated with `shipWaterForce`)|
|Cw|default = 2.5e-3|drag coefficient for whale (calculated with `whaleWaterForce`)|
|logistic||list containing `logStressCenter` and `logStressWidth` (logistic fit of an index of whale injury 0 (no injury) - 1 (fatal)), as well as `tau25`, `tau50`, `tau75` (stress that yield index values)|

In this case, we are interested in the lethality of North Atlantic right whales (NARW), but parameters on the whale body dimension and physiology can be adapted to other species.

### Quick example - 4 knots
```r
# Time interval of the simulation (1.5 sec)
t <- seq(0, 1.5, length.out = 500)

# Initial state of the model
state <- list(
  xs = -2.5,         # Ship 2.5 m away from the whale
  vs = knot2mps(4),  # 4 knot ship, use the `knot2mps()` function to convert knots to m/s
  xw = 0,            # Whale at its initial position
  vw = 0)            # Whale no moving (speed = 0)

# Parameters of the model
parms <- parameters() # Defaults

# Simulation of vessel strike
results_4kn <- strike(t, state, parms)
```
We can now plot the results of the simulation.
```r
# Plot three default graphics
par(mfcol = c(1, 3), mar = c(3.3, 3, 1, 2), mgp = c(2, 0.7, 0), cex = 0.7)
plot(results_4kn)
```
![image](https://github.com/user-attachments/assets/534662e8-f4ad-4389-9925-d39ee485115a)


These graphics show the compression of the layers of the whale. The last graph shows the lethality index over time of the collision. In this case, the strike does not reach 50% probabilty of lethality. To calculate the exact probability, we can use the `lethalityIndexFromStress()` function.
```r
# Find the maximum probability of lethality from the collision
P_leth_4kn <- max(lethalityIndexFromStress(results_4kn[["WCF"]][["stress"]]))
P_leth_4kn
```
```r
## [1] 0.295163
```
In this simulation, a collision between an adult NARW and a vessel of 45 000 kg going at a speed of 4 knots, would have a 29.5% probability of being lethal.

### Quick example - 10 knots
Let's simulate the same ship and whale but with a higher travelling speed.
```r
# Initial state of the model
state <- list(
  xs = -2.5,
  vs = knot2mps(10),  # 10 knot ship
  xw = 0,
  vw = 0)

# Simulation of vessel strike
results_10kn <- strike(t, state, parms)
```
```r
# Plot three default graphics
par(mfcol = c(1, 3), mar = c(3.3, 3, 1, 2), mgp = c(2, 0.7, 0), cex = 0.7)
plot(results_10kn)
```
![image](https://github.com/user-attachments/assets/21189493-bec9-4570-bd3e-206dd1f95fa1)

Now, we can see that the compression is more severe and that the index of lethality goes beyond 50%.
```r
# Find the maximum probability of lethality from the collision
P_leth_10kn <- max(lethalityIndexFromStress(results_10kn[["WCF"]][["stress"]]))
P_leth_10kn
```
```r
## [1] 0.6931916
```
In this simulation, a collision between an adult NARW and a vessel of 45 000 kg going at a speed of 10 knots, would have a 69.3% probability of being lethal.

### Computing probability of lethality for multiple ships going at various speed
When conducting a ship strike risk assessment, it is likely the probability of lethality needs to be calculated for multiple ships going at different speeds. For example, we have a data set of five ships of various mass and going at different speeds (these are fictional values).
| Type | Mass (kg) | Speed (kn) |
|------|--------|--------|
|Fishing|60,000|6.5|  
|Ferry|140,000|18.2| 
|Tanker|200,000|13.1|
|Cargo|250,000|9.4|
|Cruise|300,000|15.7|

We start with creating the time sequence.
```r
# Time sequence
t <- seq(0, 1.5, length.out = 500)
```

Then, we create a list to store the states of each ship based on their speed.
```r
state <- list()

for (i in 1:length(data$Type)) {
    xs = -2.5
    vs = knot2mps(vessels$Speed[i])
    xw = 0
    vw = 0
    state[[i]] <- list(xs = xs, vs = vs, xw = xw, vw = vw)
}
```
We also create a list to store the parameters of each ship based on their mass. If the impact area changes for each vessel, add the `Ly` and `Lz` arguments as well.
```r
parms <- list()

for (i in 1:length(vessels$Type)) {
    parms[[i]] <- parameters(
        ms = vessels$Mass[I])
}
```
Finally, calculate the probability of lethality for each vessel based on their unique characteristics and store the value in the data frame.
```r
# Loop to calculate the probability of lethality with speed and mass
for (i in 1:length(vessels$Type)) {
    results <- strike(t, state[[i]], parms = parms[[i]])
    vessels$Prob_Leth[i] <- max(lethalityIndexFromStress(results[["WCF"]][["stress"]]))
}
```
| Type | Mass (kg) | Speed (kn) | Probability of lethality |
|------|--------|--------|---------|
|Fishing|60,000|6.5|0.519|  
|Ferry|140,000|18.2|0.922| 
|Tanker|200,000|13.1|0.856|
|Cargo|250,000|9.4|0.749|
|Cruise|300,000|15.7|0.904|
