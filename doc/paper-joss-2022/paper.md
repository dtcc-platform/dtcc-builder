---
title: 'DTCC Builder: A mesh generator for automatic, efficient, and robust mesh generation for large-scale city modeling and simulation'
tags:
  - mesh generation
  - point cloud
  - cadastral map
  - digital twin
  - city modeling
  - city simulation
  - C++
authors:
  - name: Anders Logg
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 3"
  - name: Vasilis Naserentin
    equal-contrib: true
    affiliation: "1, 3"
  - name: Dag Wästberg
    equal-contrib: true
    affiliation: "2, 3"
affiliations:
 - name: Chalmers Univesrity of Technology
   index: 1
 - name: FIXME CIT
   index: 2
 - name: Digital Twin Cities Centre
   index: 3
date: 20 September 2022
bibliography: paper.bib
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

DTCC Builder is a mesh generator for...



DTCC Builder is implemented in C++ and very efficient...


DTCC Builder is part of DTCC Platform...


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

This work is part of the Digital Twin Cities Centre supported by
Sweden’s Innovation Agency Vinnova under Grant No.  2019-00041.

# References

# TODO

[x] The software must be open source as per the OSI definition.
[x] The software must have an obvious research application.
[x] You must be a major contributor to the software you are submitting
[x] Have a GitHub account to participate in the review process.
[x] Your paper must not focus on new research results accomplished with the software.
[x] Your paper must be hosted in a Git-based repository together with your software.
[x] Be stored in a repository that can be cloned without registration.
[x] Be stored in a repository that is browsable online without registration.
[ ] Have an issue tracker that is readable without registration.
[ ] Permit individuals to create issues/file tickets against your repository.
[ ] Indicate whether related publications exist as part of submitting to JOSS.
