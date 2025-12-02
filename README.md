# QuantizerNTDS-RD

## Description

QuantizerNTDS-RD is a project developed to illustrate the concepts presented in the paper  
**“Switched Predictor Feedbacks for Reaction-Diffusion PDEs and Globally Lipschitz Nonlinear ODE Systems Subject to Almost State and Input Quantization”** ([IMA](https://doi.org/10.1093/imamci/dnaf044)).

This project simulates:

- A globally Lipschitz nonlinear time-delay system subject to **state quantization** with a switched predictor-feedback law.
- A scalar unstable **reaction–diffusion PDE with input delay** subject to **input quantization**.

For detailed equations, assumptions, and parameter choices, please refer to the associated paper and ([technical documentation](https://tucgr-my.sharepoint.com/:b:/g/personal/fkoudohode_tuc_gr/IQB_uZ_E8hS_QKtHeHL5pfh1AVUp19CksUzti2iQVGO3v8Y?e=rnlk1t)).

## Requirements

To run this project, you will need:

- MATLAB R2023b or later.

## Installation

Follow these steps to set up the project:

1. Download the project files from [QuantizerNTDS-RD GitHub Repository](https://github.com/flo3221/quantizerNTDS-RD).
2. Extract the contents to a directory of your choice.
3. Open MATLAB and navigate to the project directory using the `cd` command:
    ```Matlab
    cd /path/to/quantizerNTDS-RD
    ```

## Usage

To use QuantizerNTDS-RD, follow these steps:

1. Open MATLAB and ensure you are in the project directory.
2. Run the main scripts:

   - For the nonlinear time-delay system with quantized predictor feedback and **dynamic** switching parameter:
     ```Matlab
     Quantized_predictor_feedback_NTDS.m
     ```

   - For the nonlinear time-delay system with quantized predictor feedback and **fixed** switching parameter:
     ```Matlab
     fixed_quantized_predictor_feedback_NTDS.m
     ```

   - For the reaction–diffusion PDE with input delay and input quantization:
     ```Matlab
     Quantization_RD_Delay.m
     ```

3. Ensure that the `private` folder is in the same directory for the nonlinear time-delay system case.  
   This folder contains routines used to solve initial–boundary value problems for first-order systems of hyperbolic partial differential equations (PDEs), following Shampine (2005).

### Functions

QuantizerNTDS-RD includes the following key functions.

#### Nonlinear time-delay system

- `hpde.m` and `setup.m`:  
  Solve the transport PDE (delay line) that appears in the ODE–PDE representation of the time-delay system.

- `mu`:  
  Implements the switching parameter $\mu(t)$ used in the dynamic quantization scheme.

- `quantizer`:  
  Implements the quantizer function for the state and input channels of the nonlinear time-delay system.

#### Reaction–diffusion PDE with input delay

- `L2norm`:  
  Computes the L²-norm  
  $\|f\|_2 = \left( \int_0^1 |f(x)|^2 \, dx \right)^{1/2}$.

- `sup_norm`:  
  Computes the sup-norm where the supremum is understood in the essential sense.


- `u0`:  
  Defines the initial condition  
    $u_0(x) = \sum_{n=1}^{3} \frac{\sqrt{2}}{n} \sin(n\pi x) + 3(x^2 - x^3)$,  
  and the corresponding constant initial condition for the delay line $v_0(x) = 5$.

- `mu_t`:  
  Implements the switching parameter $\mu(t)$ used in the quantized predictor-feedback law for the reaction–diffusion PDE.

- `quantizer_u`:  
  Implements the scalar input quantizer used for the boundary control $U(t)$.

## Examples

Refer to the following scripts for examples of how to use QuantizerNTDS-RD:

- `Quantized_predictor_feedback_NTDS.m`  
  (Nonlinear time-delay system with dynamic switching parameter.)

- `fixed_quantized_predictor_feedback_NTDS.m`  
  (Nonlinear time-delay system with fixed switching parameter.)

- `Quantization_RD_Delay.m`  
  (Reaction–diffusion PDE with input delay and input quantization.)

## Contributing

To contribute to QuantizerNTDS-RD, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or fix.
3. Make your changes and commit them.
4. Submit a pull request with a detailed description of your changes.

## License

This project is licensed under the CC BY-NC-ND license  
([`LICENSE`](https://creativecommons.org/licenses/by-nc-nd/4.0/)).

## Contact

For questions or feedback, please contact [fkoudohode@tuc.gr](mailto:fkoudohode@tuc.gr).

# Acknowledgements

Funded by the European Union (ERC, C-NORA, 101088147). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them. 

## Cite this work

If you use this code in your research, please cite:

```bibtex
@article{koudohode2025ima,
  title   = {Switched predictor feedbacks for reaction–diffusion PDEs and globally Lipschitz non-linear ODE systems subject to almost state and input quantization},
  author  = {F. Koudohode and N. Bekiaris-Liberis},
  journal = {IMA Journal of Mathematical Control and Information},
  volume  = {42},
  pages   = {1--39},
  year    = {2025}
}
