# Written Assignment 2: "Fine-tuning foundation models with MACE"

The goal of this written assignment was to "utilise machine learning interatomic potentials (MLIPs) for the fine-tuning of foundation models, specifically using the MACE architecture, [...] and to explore the advantages of fine-tuning, such as improved accuracy, reduced computational costs for training, and the effective use of smaller datasets."

This written assignment was mostly based on the following papers:

- I. Batatia et al. (2023). A foundation model for atomistic materials chemistry. arXiv preprint, 2401.00096.
- H. Kaur et al. (2025). Data-efficient fine-tuning of foundational models for first-principles quality sublimation enthalpies. Faraday Discuss., vol. 256, pg. 120&ndash;138, DOI:[10.1039/D4FD00107A](https://doi.org/10.1039/D4FD00107A).

A small collection of example datasets was also provided to students to select from.


---

In my case, I selected the dataset for defected graphene with water interfaces (for no particular reason beyond simple curiosity).


---

The project structure is as follows:

- `📁 data`: Data files to begin training from

> [!IMPORTANT]
> The actual data file, `data/input.xyz` itself, is not in this repository due to large file size!

- `📁 finetune`: Fine-tuned model, derived from MACE-MP-0

- `📁 scratch`: Model trained from scratch, with the same hyperparameters as the finetuned model

- `📁 eval`: Series of tests and evaluations, to compare between the two models

