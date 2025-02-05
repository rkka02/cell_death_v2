# Cell death project (2024.06. ~ ???)
> The ikigai for CDD, since its first conception, has been the advancement of quality science.
> This is something that gives CDD a sense of purpose, a reason for living, and a commitment to serving both authors and readers.
> (Nature Cell Death & Differentiation)[^1]


## What:
Develop an AI model which can determine whether a cell is dead or alive. Accuracy goal is over 95%. IF this is cleared, then we can further specify cell death types : Necrosis, Necroptosis, Apoptosis, and Pyroptosis.

## Why:
Bio companies these days do millions of experiments parallely in a day. Human can't check each result manually, so automation of it will be meaningful.

## How:
Construct a deep learning model which predicts whether a cell is dead or alive. Input data is a time snapshot of a 3D FOV(cells). We won't use time-lapsed. It's another story.

## Problems


## Study
After reading a review paper on QPI technique[^2], I could gather fundamentally novel parameters that QPI can provide. They are 
- Cell (dry) mass density
- 3D morphological parameters : cell volume, surface area, sphericity
- cell membrane deformability

[^1]: Melino, G., Knight, R.A., Mak, T.W. et al. The birth of death, 30 years ago. Cell Death Differ 31, 379–386 (2024). https://doi.org/10.1038/s41418-024-01276-8
[^2]: Park, Y., Depeursinge, C. & Popescu, G. Quantitative phase imaging in biomedicine. Nature Photon 12, 578–589 (2018). https://doi.org/10.1038/s41566-018-0253-x
