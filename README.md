# FloodRiskManagement_model_Rhine: A python-based flood risk assessment model for the lower Rhine River.
The model allows estimating flood risk accounting for hydrological and dike breaching uncertainties for the Lower Rhine River (see in Figure 1), involving Germany and The Netherlands.

|<img src="figs/Rhine.png" width="600"/>|
| ------------- |
| *Figure 1: Study area. Six areas are identified: four Dutch areas (area 0, in red; area 1, in blue; area 2, in orange; 114 area 3, in black) and two German areas (area 4, in green; area 5, in purple). The administrative country 115 border is depicted in grey.*          |

The modeling steps are summarized in Figure 2. Flood routing is modelled through a Muskingum scheme which has been benchmarked on results from a SOBEK model calibrated on the case study. Three flood risk reduction measures are possible: embankment heightening, making room for the river and flow diversion at the bifurcation points. Embankments can be raised up to 1 meter, with steps of 10 centimeters. 156 Room for the River projects are available along the Dutch Rhine. The implementation of each project is established through a binary decision variable, which is 1 when it is implemented and 0 otherwise. There are two bifurcation points, and the default flow distributions are the ones provided by the SOBEK model. At each point, a distribution change of plus/minus 30 % of the default distribution can be implemented. 
Outputs in terms of expected annual damage relate to each of the dike ring area; investment costs relate to the different interventions being adopted.

|<img src="figs/model_scheme.png" width="600"/>|
| ------------- |
| *Figure 2: The modeling scheme*  |
