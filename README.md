# tess_gaia_luminosity
Compute TESS luminosity using colors measured with Gaia

## Example usage


````python
from TESSgaiaLum import getTESSlum

TIC = 267166222
getTESSlum(TIC)

list_of_tics = [267166222, 323433390]
getTESSlum(list_of_tics)
````
