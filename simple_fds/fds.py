import os
import sys
import argparse
from copy import deepcopy

import numpy as np

from .timedilation import TimeDilation, TimeDilationExact
from .segment import run_segment, trapez_mean
from .lsstan import LssTangent#, tangent_initial_condition
from .timeseries import windowed_mean
# ---------------------------------------------------------------------------- #

def tangent_initial_condition(subspace_dimension,total_dimension):
    W = np.random.rand(subspace_dimension,total_dimension)
    W = np.linalg.qr(W.transpose())[0].transpose()
    w = np.zeros(total_dimension)
    return W, w

def lss_gradient(lss, G_lss, g_lss, G_dil, g_dil, segment_range=None):
    if segment_range is not None:
        lss = deepcopy(lss)
        if isinstance(segment_range, int):
            s = slice(segment_range)
        else:
            s = slice(*segment_range)
        lss.bs = lss.bs[s]
        lss.Rs = lss.Rs[s]
        G_lss  = G_lss [s]
        g_lss  = g_lss [s]
        J      = J     [s]
        G_dil  = G_dil [s]
        g_dil  = g_dil [s]
    alpha = lss.solve()
    grad_lss = (alpha[:,:,np.newaxis] * np.array(G_lss)).sum(1) + np.array(g_lss)
    J = np.array(J)
    dJ = trapez_mean(J.mean(0), 0) - J[:,-1]
    steps_per_segment = J.shape[1]
    dil = ((alpha * G_dil).sum(1) + g_dil) / steps_per_segment
    grad_dil = dil[:,np.newaxis] * dJ
    return windowed_mean(grad_lss) + windowed_mean(grad_dil)

class RunWrapper:
    def __init__(self, run):
        self.run = run

    def variable_args(self, u0, parameter, steps):
        return self.run(u0, parameter, steps)
    def __call__(self, u0, parameter, steps):
        u1, J = self.variable_args(
                    u0, parameter, steps)
            return u1, np.array(J).reshape([steps, -1])
      
def continue_shadowing(
        run, parameter,
        num_segments, steps_per_segment, epsilon=1E-6,V, v, u0, lss):
    """
    """
    

    run = RunWrapper(run)
    u0, V, v, lss, G_lss, g_lss, J_hist, G_dil, g_dil = checkpoint


    i = lss.K_segments()
    time_dil = TimeDilation(run, u0, parameter)

    V = time_dil.project(V)
    v = time_dil.project(v)

    u0, V, v, J0, G, g = run_segment(
            run, u0, V, v, parameter, i, steps_per_segment,
            epsilon)

    J_hist.append(J0)
    G_lss.append(G)
    g_lss.append(g)

    for i in range(lss.K_segments() + 1, num_segments + 1):

        # time dilation contribution
        run_id = 'time_dilation_{0:02d}'.format(i)
        if run_ddt is not None:
            time_dil = TimeDilationExact(run_ddt, u0, parameter)
        else:
            time_dil = TimeDilation(run, u0, parameter, run_id,
                                    simultaneous_runs, interprocess)
        G_dil.append(time_dil.contribution(V))
        g_dil.append(time_dil.contribution(v))

        V = time_dil.project(V)
        v = time_dil.project(v)

        V, v = lss.checkpoint(V, v)
        # extra outputs to compute
        compute_outputs = [lss.Rs[-1], lss.bs[-1], G_dil[-1], g_dil[-1]]

        # run all segments
        if i < num_segments:
            u0, V, v, J0, G, g = run_segment(
                    run, u0, V, v, parameter, i, steps_per_segment,
                    epsilon, simultaneous_runs, interprocess, get_host_dir=get_host_dir,
                    compute_outputs=compute_outputs, spawn_compute_job=spawn_compute_job)
        else:
            run_compute(compute_outputs, spawn_compute_job=spawn_compute_job, interprocess=interprocess)

        for output in [lss.Rs, lss.bs, G_dil, g_dil]:
            output[-1] = output[-1].field

        checkpoint = Checkpoint(
                u0, V, v, lss, G_lss, g_lss, J_hist, G_dil, g_dil)
        print(lss_gradient(checkpoint))
        sys.stdout.flush()

        if checkpoint_path and (i) % checkpoint_interval == 0:
            save_checkpoint(checkpoint_path, checkpoint)

        if i < num_segments:
            J_hist.append(J0)
            G_lss.append(G)
            g_lss.append(g)

    if return_checkpoint:
        return checkpoint
    else:
        G = lss_gradient(checkpoint)
        return np.array(J_hist).mean((0,1)), G

def shadowing(
        run, u0, parameter, subspace_dimension, num_segments,
        steps_per_segment, runup_steps, epsilon=1E-6):
    '''
    run: a function in the form
         u1, J = run(u0, parameter, steps)

         inputs  - u0:           init solution, a flat numpy array of doubles.
                   parameter:    design parameter, a single number.
                   steps:        number of time steps, an int.
         outputs - u1:           final solution, a flat numpy array of doubles,
                                 must be of the same size as u0.
                   J:            quantities of interest, a numpy array of shape
                                 (steps, n_qoi), where n_qoi is an arbitrary
                                 but consistent number, # quantities of interest.
    '''
    

    run = RunWrapper(run)
  
    if runup_steps > 0:
        u0, _ = run(u0, parameter, runup_steps)

    total_dimension = len(u0) 
    V, v = tangent_initial_condition(subspace_dimension,total_dimension)
    lss = LssTangent()
    return continue_shadowing(
            run, parameter,
            num_segments, steps_per_segment, epsilon, V, v, u0, lss)
