# Liquid State Machine for Speech Classification

This repository contains the MATLAB implementation of a Liquid State Machine (LSM) for speech classification, demonstrating the use of an approximate state-space model for predictive performance, as described in the papers listed in the citation section.

## Citing

If you find this repository useful in your research, please consider citing the following papers:

```bibtex
@inproceedings{gorad2019predicting,
  title={Predicting performance using approximate state space model for liquid state machines},
  author={Gorad, Ajinkya and Saraswat, Vivek and Ganguly, Udayan},
  booktitle={2019 International Joint Conference on Neural Networks (IJCNN)},
  pages={1--8},
  year={2019},
  organization={IEEE}
}

@inproceedings{saraswat2021hardware,
  title={Hardware-friendly synaptic orders and timescales in liquid state machines for speech classification},
  author={Saraswat, Vivek and Gorad, Ajinkya and Naik, Anand and Patil, Aakash and Ganguly, Udayan},
  booktitle={2021 International Joint Conference on Neural Networks (IJCNN)},
  pages={1--8},
  year={2021},
  organization={IEEE}
}
```
## Usage

Set the simulation parameters in ExamplePrescript.m.
Run the main.m script to start the LSM simulation, which in turn calls SpokenDigitsLSM.m.
### Prerequisites
MATLAB R2021a or later
Auditory Toolbox for MATLAB (included in the repository)
### Installation
Clone the repository and run the main.m script in MATLAB.

### Output
The simulation outputs results in a .mat file along with log files. 
