# Optimization Framework for Selection and Scaling of Earthquake Records

## Workflow
![image](https://user-images.githubusercontent.com/61163577/212420429-8876eeb2-5319-48f0-9d47-522cf824a0dc.png)

## Example
```  
import scalepy

target = scalepy.targetSpectrum(0.8, 0.4, 'ZD')
keys = scalepy.recordSelection('4 9', '180 360', '0 20', 'Strike - Slip', '0 50', '0 50', '0 50', target, period=1)
scalepy.amplitudeScaling(keys, target, 1, spectral_ordinate='srss')
```

https://dogukankaratas-scalepy-guiscalepygui-5jmoeu.streamlit.app/
 
## License
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
