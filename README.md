# DCAXL

The script provided herein (sbmbuild.py) can be used to generate structure-based models (SBMs) from interaction data obtained from coevolution analysis and chemical cross-linking/SM and run protein folding predictions.

The detailed description of application of these strategies to model building is described in:

R. N. dos Santos, A. J. R. Ferrari, H. C. R. de Jesus, F. C. Gozzo, F. Morcos, L. Martínez, Enhancing protein fold determination by exploring the complementary information of chemical cross-linking and coevolutionary signals. Bioinformatics, 2018
https://doi.org/10.1093/bioinformatics/bty074

The following input data should be provided:

    (1) Primary sequence of target protein.
  
    (2) Predicted type of secondary structure (SS) for each sequence region. Any SS prediction server can be used (e.g Jpred,       
      PSIPRED, PP, etc). However, SS tags should be represented as follows: "H" for alpha-helices, "E" for beta-strands and "C" 
      or "-" for undetermined (coiled-coil) structures. 
      
    (3) List of high-correlated residue couplings estimated from coevolutionary analysis. Despite any available method can be 
      used, we highly recommend Direct-Coupling Analysis approach (freely available at http://dca.rice.edu/portal/dca/).
      
    (4) List of residue pairs identified by chemical cross-linking/mass spectrometry experiments. The most common chemical 
      linkers are supported: 
            (A) Direct reaction (zero-length) between residue sidechains in the presence of carbodiimides. 
                Possible pairs: DK, EK, DS, ES.
            (B) Disuccinimidyl suberate (DSS). Possible pairs: KK, KS, SS.
            (C) 1,6-Hexanediamine. Possible pairs: EE, ED, DD
	
Data for (1) and (2) should be compiled in a text file with sequence and predicted SS located in first and second lines, respectively (as shown above).

        ...EDDMEVVGKGTTGQDAVDFVKKRQP...
        ...CCCCEEEEEECCHHHHHHHHHHCCC...

Finally, an initial 3D model and a topology can be generated using the following command:
  
        python sbmbuild.py sequence_SS coevolution_list crosslink_list

Input/output examples are found herein (.TAR.GZ files). 
Fold predictions can be carried in GROMACS/gaussian simulations (http://smog-server.org/extension/gromacs-4.5.4_sbm1.0.tar.gz) using the generated .gro and .top files. A example GROMACS parameter file for an annealing protocol based in this proposed methodology is also provided and can be used for any system of study (sbm_calpha_SA_200to1_20ns.mdp).

Further assistance can be found at: rikchicfb@gmail.com
