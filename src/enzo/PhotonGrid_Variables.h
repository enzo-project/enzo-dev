#define MAX_SOURCES 1000
// Photons
int    NumberOfPhotonPackages;

//  start pointer for linked list of packages
PhotonPackageEntry *PhotonPackages;

// linked list of packages with its work already finished
PhotonPackageEntry *FinishedPhotonPackages;  

// linked list of packages that are "paused", waiting to be merged
PhotonPackageEntry *PausedPhotonPackages;

// Sources
//int    NumberOfRadiationSources;
grid  **SubgridMarker; // pointers to first array elements of subgrids for each cell
int HasRadiation;
