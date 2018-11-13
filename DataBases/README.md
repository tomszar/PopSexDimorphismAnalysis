# Description of databases

- Covariates.csv: contains demographic and anthropometric information from participants, including ID, Sex, Age, Height, and Weight
- FacePCA: contains databases obtained from the geometric morphometrics pipeline
    + cs.csv: centroid size of each face
    + eigenvalues.csv: eigenvalues of the 100 principal components extracted from the landmark data
    + facets.csv: contains the triangulation information to generate the 3D surface from 3D landmarks
    + landmark_ids.csv: contains the ID information
    + scores.csv: PCA scores for each participant
- Genotypes: contains databases obtained from the genotype pipeline
    + geno.eigenval: eigenvalues from PCA on genotype data
    + geno.eigenvec: eigenvectors from PCA on genotype data
    + geno.fam: contains the ID information of genotype files
    + RefPops.csv: contains the ID, pop, and super_pop information from reference genotype data
