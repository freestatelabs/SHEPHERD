# SHEPHERD
Open-source mechanical FEA that *just works*.

* High-performance structural and thermal calculations for the design engineer
* Easily deploy on Windows, Linux, and MacOS - all required libraries are included in this repository
* Highly parallelized operations, use every processor core available
* Focused feature set; written with eye towards simplicity and maintainability
* Thoroughly documented code and theory manual with validation examples
* Interoperability with Calculix input files
* Static, modal, transient, & nonlinear analysis types all supported

# Feature List 

## Tier 1 
* Structural analysis
  * Linear static
  * Nonlinear static
  * Features
    * Large deflection
    * Contact
      * Sliding
      * Bonded
    * Loads
    * Boundary Conditions
  * Buckling
  * Modal
  * Sine vibe
  * Random vibe
  * Transient
    * Mode superposition
    * Full transient
* Linear static thermal analysis 
* Element library (first and second order)
  * 1D elements
    * Beam
    * Spring (1-3 DOF)
  * 2D elements
    * Tri3
    * Quad4
  * 3D elements
    * Tet4
    * Ted8
    * Py5
    * Pyr13
    * Wed6
    * Wed15
    * Hex8
    * Hex20
* Solvers
  * Time
    * Static
    * Transient
  * Frequency
    * Eigenmodes
    * Mode superposition
  * Linear
  * Nonlinear