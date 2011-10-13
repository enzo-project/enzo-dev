/******************************
/
/ List of Photons to be moved
/
*****************************/ 

struct ListOfPhotonsToMove {
  PhotonPackageEntry  *PhotonPackage;
  ListOfPhotonsToMove *NextPackageToMove;
  grid *FromGrid;
  grid *ToGrid;  
  int ToLevel;
  int ToProcessor;
  int ToGridNum;
  char PausedPhoton;
};
