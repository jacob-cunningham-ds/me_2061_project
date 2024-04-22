# Field Reconstruction for Unseen Data

In this case, we assume that we have temperature sensors at the top surface of the cylinder ($0 \leq \theta \leq 180$) for a new (unknown) heat source. The objective is to use the measurements of these sensors to reconstruct the temperature field everywhere in the domain.

You need to first generate the true temperature field for generating synthetic data for sensors as well for verification purposes. To this end, run the FOM for $\theta = 15$. Use the `ExtractCylinderT` function to extract the temperature at the top surface of the cylinder. Assume that these are the sensor measurements.

Use these sensor measurements along with the POD modes to reconstruct the temperature at all points. Plot the reconstructed temperature field and true temperature field. Also, plot the contours of error by subtracting the reconstructed field from the true field.

MATLAB code was used to reconstruct the temperature field everywhere in the domain, run the FOM for $\theta = 15$, use `ExtractCylidnerT` function to extract the temperature at the top surface of the cylinder, plot the reconstructed temperature field and true temperature field as shown in Figure BLANK, and plot the contours of error as shown in {numref}`field_compare`.

## Temperature Field Comparison

```{figure} img/field_compare.png
---
name: field_compare
---
Comparison of the temperature fields at $\theta = 15^{o}$
```