Intrinsic Signal Imaging to compute ocular dominance


Images are acquired using the ContImage program that runs a Dalsa 1M30 frame transfer CCD camera under Windows.  The version of ContImage used was V1.01 from March 20, 2006.  Code was compiled using Microsoft Visual C++.  ContImage was written by Valery Kalatsky.  It depends on the PGPLOT graphics subroutine library https://sites.astro.caltech.edu/~tjp/pgplot/  Image acquisition is described in 
Kalatsky VA, Stryker MP. New paradigm for optical imaging: temporally encoded maps of intrinsic signal. Neuron. 2003 May 22;38(4):529-45. doi: 10.1016/s0896-6273(03)00286-1. PMID: 12765606,
avalable at https://www.sciencedirect.com/science/article/pii/S0896627303002861?via%3Dihub


Images were analyzed using the Iman suite written by Valery Kalatsky, also described above.  Documents are in ISI_ImAnManual.pdf and ISI_cheatsheet.pdf


Visual stimuli for the mice were generated on an Ubuntu linux computer with duplicate screens (one for the operator and one for the mouse)  using the standard Matlab Psychtoolbox libraries by the PsychStimController script running with  the ContStim UDP synchronization.  PsychStimController was written by Cristopher Niell and Michael Stryker.


Ocular dominance measurements were computed in Matlab using the scripts in JC-ImanMaps written mostly by Jianhua Cang, as described in
Cang J, Kalatsky VA, Löwel S, Stryker MP. Optical imaging of the intrinsic signal as a measure of cortical plasticity in the mouse. Vis Neurosci. 2005 Sep-Oct;22(5):685-91. doi: 10.1017/S0952523805225178. PMID: 16332279; PMCID: PMC2553096.