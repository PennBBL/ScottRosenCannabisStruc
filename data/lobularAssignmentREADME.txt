ROI information for Scott et al. Scott et al. exploring relationship of cannabis consumption with structural neuroimaging data. Images used for the fusion process, to determine gyral and sulcal landmarks, can be found by following the link provided: https://my.vanderbilt.edu/masi/about-us/resources-data/

If using these labels please cite: Yuankai Huo, Andrew J. Asman, Andrew J. Plassard, Bennett A. Landman. "Simultaneous total intracranial volume and posterior fossa volume estimation using Non-local Spatial STAPLE label fusion." Human Brain Mapping. 2016

Furthermore, structural processing providing the regions used in Scott et al. can be performed using the xcpEngine. The xcpEngine can be found here: https://github.com/PennBBL/xcpEngine

If using the xcpEngine please cite: Ciric R, Rosen AFG, Erus G, Cieslak M, Adebimpe A, Cook PA, Bassett DS, Davatzikos C, Wolf DH & Satterthwaite TD. "Mitigating head motion artifact in functional connectivity MRI." Nature Protocols 1 (2018). doi:10.1038/s41596-018-0065-y

The lobularAssignment.csv contains information about regions included in all analyses.

Column definitions can be found below:

ROI_readable: Full name of region
ROI: abbreviated name of region
Lobular_Assignment: Details which ROI's were included in the lobes. To create the volume lobes, regions were summed, cortical thickness and GMD take the volume weighted mean across all corresponding regions
Lateral: Details if there is a lateral component to the identified region (Y=Yes, N=No)
Vol: Details if this region was included when comparing volumes across groups (I=Included, E=Excluded)
CT: Details if this region was included when comparing cortical thickness across groups (I=Included, E=Excluded)
GMD: Details if this region was included when comparing gray matter density across groups (I=Included, E=Excluded)
