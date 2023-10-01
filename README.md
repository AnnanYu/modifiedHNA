# Modified HNA Algorithm

This is the code repository accompanied the manuscript titled "Leveraging the Hankel norm approximation and data-driven algorithms in reduced-order modeling." It implements the two-stage algorithm described in the paper by combing the regularized block-AAA algorithm and the modified HNA algorithm. The block-AAA part was heavily adapted from the code in the [block_AAA](https://github.com/nla-group/block_aaa) repository. The benchmark examples were taken from the [SLICOT](http://www.slicot.org/) benchmark collection. The code was implemented using MATLAB 2023b. In addition, the [Control System Toolbox](https://www.mathworks.com/products/control.html) is needed.

## Repository structure
Every MATLAB file starting with the word `Driver` is used to reproduce an experiment in the manuscript. The functions `block_aaa` and `block_aaa_disk` contain implementation of the block-AAA part. The function `modified_HNA` contains the implementation of the modified HNA part. Every other MATLAB file is a helper function.

## Executing examples in the manuscript
Every experiment can be reproduced using a MATLAB script.

* To reproduce the experiment in section 4.4, run `Driver_stability.m`.
* To reproduce the experiment in section 5.4, run `Driver_toy_hilbert.m`.
* To reproduce the experiment in section 7.1, run `Driver_modified_HNA.m`.
* To reproduce the experiment in section 7.2, run `Driver_eady.m`.
* To reproduce the experiment in section 7.3, run `Driver_CD_player.m`.

The scripts `Driver_eady.m` and `Driver_CD_player.m` employs the block-AAA algorithm with a regularizer. To disable the regularizer, set the regularizing parameter to zero.
