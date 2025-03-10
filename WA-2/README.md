# Written Assignment 2: "Fine-tuning foundation models with MACE"

The goal of this written assignment was to "utilise machine learning interatomic potentials (MLIPs) for the fine-tuning of foundation models, specifically using the MACE architecture, [...] and to explore the advantages of fine-tuning, such as improved accuracy, reduced computational costs for training, and the effective use of smaller datasets."

This written assignment was mostly based on the following papers:

- I. Batatia et al. (2023). A foundation model for atomistic materials chemistry. arXiv preprint, 2401.00096.
- H. Kaur et al. (2025). Data-efficient fine-tuning of foundational models for first-principles quality sublimation enthalpies. Faraday Discuss., vol. 256, pg. 120&ndash;138, DOI:[10.1039/D4FD00107A](https://doi.org/10.1039/D4FD00107A).

A small collection of example datasets was also provided to students to select from.


---

In my case, I selected the dataset for (up to) one molecule of NaCl dissolved in water, for no particular reason beyond simple curiosity.


---

The project structure is as follows:

- `📁 venv`: Python virtual environment needed to run the MACE models

> [!IMPORTANT]
> This repository does not contain the Python environment files themselves, due to large filesizes! Run the script `venv/setup.sh` to download and install the environment, and `source venv/bin/activate` to activate the environment, before running any script.

- `📁 data`: Dataset to train on

- `📁 finetune_naive`: Fine-tuned model, derived from MACE-MP-0b, using naive transfer learning

- `📁 finetune_multi`: Fine-tuned model, derived from MACE-MP-0b, using multihead replay to prevent catastrophic forgetting

- `📁 scratch`: Model trained from scratch, with the same hyperparameters as the finetuned model

- `📁 eval`: Series of tests and evaluations, to compare between the two models

