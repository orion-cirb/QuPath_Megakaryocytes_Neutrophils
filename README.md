# QuPath_Megakaryocytes_Neutrophils

* **Developed for:** Cécile
* **Team:** De Thé
* **Date:** June 2023
* **Software:** QuPath

### Images description

2D images of mouse bone marrow taken with the Axioscan

3 channels: 
  1. *DAPI:* nuclei
  2. *EGFP:* megakaryocytes
  3. *Cy5:* neutrophils (not mandatory)
  

### Plugin description

* In each annotation, detect tissue with a pixel classifier
* In tissue, detect megakaryocytes with Cellpose
* If Cy5 channel exits, detect neutrophils with Cellpose
* Only keep neutrophils in megakaryocytes

### Dependencies

* **QuPath pixel classifier** named *Tissue.json*
* **Cellpose** QuPath extension and conda environment + *cyto2* model

### Version history

Version 1 released on June 14, 2023.
