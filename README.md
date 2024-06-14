# The Mixture Adaptive Design (MAD): An Experimental Design for Anytime-Valid Causal Inference on Multi-Armed Bandits
 The Mixture Adaptive Design (MAD) is a new experimental design for multi-armed bandit (MAB) algorithms that enables anytime valid inference on the average treatment effect for any MAB algorithm. Intuitively, the MAB "mixes" any bandit algorithm with a Bernoulli design such that at each time step, the probability that a customer is assigned via the Bernoulli design is controlled by a user-specified deterministic sequence that can converge to 0. The sequence enables managers to directly and interpretably control the trade-off between regret minimization and inferential precision. Under mild conditions on the rate the sequence converges to 0, we provide a confidence sequence that is asymptotically anytime valid and demonstrate that the MAD is guaranteed to have a finite stopping time in the presence of a true non-zero ATE. Hence, the MAD allows managers to stop experiments early when a significant ATE is detected while ensuring valid inference, enhancing both the efficiency and reliability of adaptive experiments. Empirically, we demonstrate that the MAD achieves finite-sample anytime-validity while accurately and precisely estimating the ATE, all without incurring significant losses in reward compared to standard bandit designs.

 For details about the methodology, check out our paper:
https://arxiv.org/abs/2311.05794

This repository contains the simulations.R code file, which generates the data for reproducing all simulation results in the paper. 

