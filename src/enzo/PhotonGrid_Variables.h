#define MAX_SOURCES 1000
// Photons
int    NumberOfPhotonPackages;
int    NumberOfRenderingPackages;

//  start pointer for linked list of packages
PhotonPackageEntry *PhotonPackages;

// linked list of packages with its work already finished
PhotonPackageEntry *FinishedPhotonPackages;  

// linked list of packages that are "paused", waiting to be merged
PhotonPackageEntry *PausedPhotonPackages;

// linked list of packages used for projections or volume renderings
PhotonPackageEntry *RenderingPhotonPackages;

// Sources
//int    NumberOfRadiationSources;
grid  **SubgridMarker; // pointers to first array elements of subgrids for each cell
int HasRadiation;
