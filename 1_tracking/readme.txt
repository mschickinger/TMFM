This folder contains code for the localization of fluorescently-labeled objects over the course of TIRF microscopy recordings (=position tracking). More specifically, I used it to track the positions of the reporter units in TMFM.

The main script is the ’TWO_COLOR_TRACKER_v13c.m’. All the other functions are called from within that script. It should be noted that this script is written for execution on a computing cluster with a large file-storage system (nfs). The frame-by-frame evaluation of the recordings with e.g. 45000 frames in two spectral channels was performed in parallel with up to 64 workers. Given enough time, same procedure would of course also be feasible on a desktop machine. The code would then have to be modified accordingly.

The foundation for the whole analysis process is the construction of a so-called “movie”-object for every dataset. I would like to acknowledge the major contribution of my colleague Jonas Funke to the ‘movie.m’ class definition file. In fact, it was he who took the initiative to “go object-oriented” in the first place.

For estimating the centroid positions of the objects, I adapted the VWCM algorithm originally published in: A. J. Berglund, M. D. McMahon, J. J. McClelland, and J. A. Liddle, “Fast, bias-free algorithm for tracking single particles with variable size and shape,” Optics Express, vol. 16, no. 18, pp. 14064–14075, 2008.

- - - - -
Matthias Schickinger
October 2020