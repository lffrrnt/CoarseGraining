# CoarseGraining

Code accompanying the preprint:

> "Scaling and tuning to criticality in resting-state human magnetoencephalography"
> [Preprint link](https://arxiv.org/pdf/2602.17820)
>
> ## Abstract
> From 1/f noise to neuronal avalanches, evidence of scaling in brain activity has been increasingly linked to tuning to or near criticality. The concept of scaling is intimately related
>  to the renormalization group (RG), in essence providing coarse-grained, simplified descriptions that generalize to classes of diverse physical systems. Following the RG idea,
> scaling laws have been reported in populations of spiking neurons at microscopic scales. Whether similar scaling principles govern large-scale neural activity in the human brain
>  and how they relate to underlying neural physiology remain unresolved. Here, we analyze large-scale electrophysiological recordings (MEG) of human resting-state brain activity
> and apply a RG-inspired coarse-graining approach to track collective neural dynamics across spatial scales. We find that multiple observables exhibit robust scale-invariant behavior under
> coarse-graining: activity variance and correlations grow according to power laws, covariance eigenspectra follow a characteristic scaling relation, and neuronal avalanche statistics
> remain invariant. Using an analytically tractable neural network model, we show that the observed scaling signaturesarise when the system operates slightly below criticality, and that
> the scaling exponents depend on the excitation–inhibition balance. These findings demonstrate that RG-inspired scaling analysis can uncover signatures of critical dynamics in
> non-invasive human electrophysiology and suggest a principled route toward estimating excitation–inhibition balance from large-scale brain recordings.

This repository contains Python scripts implementing the coarse-graining procedure introduced by Meshulam et al. in 
*"Coarse Graining, Fixed Points, and Scaling in a Large Population of Neurons"*. These scripts were used for the coarse-graining analysis described in our manuscript/preprint.

## Files

- `cg_funcs.py` performs iterative coarse-graining steps and computes cluster-level covariance eigenspectra.
- `analyse_quantities.py` include functions for computing silence probability, cluster variance, and autocorrelation properties of coarse-grained variables across renormalization steps.

## Requirements

Python 3 with standard scientific libraries:

- numpy
- scipy
- matplotlib

## Citation

If you use this code, please use:

@article{https://arxiv.org/pdf/2602.17820,
  title={Scaling and tuning to criticality in resting-state human magnetoencephalography},
  author={Topal, Irem and Poggialini, Anna and Maschio, Marco Dal and De Martino, Daniele and Shriki, Oren and Lombardi, Fabrizio},
  journal={arXiv preprint arXiv:2602.17820},
  year={2026}
}

## License

This project is licensed under the MIT License.
