# Optimization Framework for Selection and Scaling of Earthquake Records

## Workflow
![flowchart](https://user-images.githubusercontent.com/61163577/159950981-13642764-8c2b-4952-9168-aef0bc762b33.png)

## Example
```  
import scaleopt

target = scaleopt.targetSpectrum(0.8, 0.4, 'ZD')
keys = scaleopt.recordSelection('4 9', '180 360', '0 20', 'Strike - Slip', '0 50', '0 50', '0 50', target, period=1)
scaleopt.amplitudeScaling(keys, target, 1, spectral_ordinate='srss')
```
 

## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)