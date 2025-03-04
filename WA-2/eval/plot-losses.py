#!/usr/bin/env python

import json
import numpy as np
import matplotlib.pyplot as plt

#==================================================================================================
# Script parameters

INPUT_FILES = {'scratch': '../scratch/results/model_from_scratch_run-123_train.txt',
               'finetune_naive': '../finetune_naive/results/model_finetuned_naive_run-123_train.txt',
               'finetune_multihead': '../finetune_multihead/results/model_finetuned_naive_run-123_train.txt'} # N.B. there was mistake in the model name for this one

#==================================================================================================
# Do not modify below this line

fig, axes = plt.subplots(2, 2)

for label, filename in INPUT_FILES.items():

    eval_loss_over_epoch = list()
    opt_loss_over_epoch = list()
    rmse_e_over_epoch = list()
    rmse_f_over_epoch = list()
    opt_data = list()
    current_epoch = 0
    with open(filename, 'r') as f:
        for line in f:
            fields = json.loads(line)
            if 'epoch' in fields and fields['epoch'] is not None:
                if fields['epoch'] == current_epoch:
                    if fields['mode'] == 'opt':
                        opt_data.append(fields['loss'])
                    elif fields['mode'] == 'eval':
                        if current_epoch == len(eval_loss_over_epoch):
                            eval_loss_over_epoch.append(fields['loss'])
                            rmse_e_over_epoch.append(fields['rmse_e_per_atom'])
                            rmse_f_over_epoch.append(fields['rmse_f'])
                        else:
                            eval_loss_over_epoch[-1] = (eval_loss_over_epoch[-1] + fields['loss']) / 2
                            rmse_e_over_epoch[-1] = (eval_loss_over_epoch[-1] + fields['rmse_e_per_atom']) / 2
                            rmse_f_over_epoch[-1] = (eval_loss_over_epoch[-1] + fields['rmse_f']) / 2
                else:
                    opt_loss_over_epoch.append(np.mean(opt_data))
                    opt_data = list()
                    current_epoch = fields['epoch']
                    opt_data.append(fields['loss'])

    axes[0][0].plot(eval_loss_over_epoch, label=label)
    axes[1][0].plot(opt_loss_over_epoch, label=label)
    axes[0][1].plot(rmse_e_over_epoch, label=label)
    axes[1][1].plot(rmse_f_over_epoch, label=label)

axes[0][0].set_title('Loss (during evaluation) over epoch number')
axes[0][0].set_xlabel('Epoch no.')
axes[0][0].set_ylabel('Loss')
axes[0][0].set_yscale('log')
axes[0][0].legend()

axes[1][0].set_title('Mean loss (during optimization) over epoch number')
axes[1][0].set_xlabel('Epoch no.')
axes[1][0].set_ylabel('Mean loss')
axes[1][0].set_yscale('log')
axes[1][0].legend()

axes[0][1].set_title('RMSE energy per atom over epoch number')
axes[0][1].set_xlabel('Epoch no.')
axes[0][1].set_ylabel(r'RMSE energy per atom ($eV$)')
axes[0][1].set_yscale('log')
axes[0][1].legend()

axes[1][1].set_title('RMSE force over epoch number')
axes[1][1].set_xlabel('Epoch no.')
axes[1][1].set_ylabel(r'RMSE force ($eV/\AA$)')
axes[1][1].set_yscale('log')
axes[1][1].legend()

plt.show()


