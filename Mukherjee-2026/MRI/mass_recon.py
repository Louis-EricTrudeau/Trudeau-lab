import os
import shutil
import sys
import getopt
from bruker2nifti.converter import Bruker2Nifti

"""

Ensemble reconstruction of several experiments that are contained in a single directory.
Requires the Bruker2Nifti package

Usage:
$ conda activate <your_env>
$ pip install Bruker2Nifti
$ python mass_recon.py -i <path/to/dir -o <reconstruction_folder_name>

Quantitative Microstructure Imaging Lab
Written by VG 2022


"""


def rename_scan_dirs(input_study_dir, recond_folder_name):
    """ assign reconstructed scan directories with meaningful names based on acquisition method """
    scan_dirs = sorted(os.listdir(os.path.join(input_study_dir, recond_folder_name)))
    study_dirs_abs = os.path.join(input_study_dir, recond_folder_name)
    # cd to absolute study directory and create temporary folder
    os.chdir(study_dirs_abs)
    os.mkdir('temp')

    # loop over all scans within current study
    for curr_scan_idx in range(scan_dirs.__len__()):
        scan_dir_abs = os.path.join(study_dirs_abs, scan_dirs[curr_scan_idx])
        # open directory, read in acquisition method file
        os.chdir(scan_dir_abs)
        method_file = open('acquisition_method.txt', mode='r')
        method_txt = method_file.read()
        method_file.close()
        print(method_txt)
        for file in os.listdir():
            if file.startswith(recond_folder_name):
                os.rename(file, file.replace(recond_folder_name, method_txt))

        # move files out of directory, rename it, then move the files back in
        for file in os.listdir():
            shutil.move(os.path.join(os.getcwd(), file), os.path.join(study_dirs_abs, 'temp'))

        os.rename(os.path.join(study_dirs_abs, os.path.basename(os.getcwd())),
                  os.path.join(study_dirs_abs,  os.path.basename(os.getcwd()).replace(recond_folder_name, method_txt)))

        for file in os.listdir(os.path.join(study_dirs_abs, 'temp')):
            shutil.move(os.path.join(os.path.join(study_dirs_abs, 'temp'), file),  os.getcwd())

    # delete temporary directory
    os.rmdir(os.path.join(study_dirs_abs, 'temp'))
    return 0


def reconstruct_study(input_study_dir, recond_folder_name):
    """ reconstruct all scans within a specified study using the Bruker2Nifti package """
    # specify directories
    pfo_study_in = input_study_dir
    pfo_study_out = pfo_study_in
    # instantiate a Bruker Reconstructor object and its parameters
    bru = Bruker2Nifti(pfo_study_in, pfo_study_out, study_name=recond_folder_name)
    bru.verbose = False
    bru.correct_slope = True
    bru.get_acqp = True
    bru.get_method = True
    bru.get_reco = False
    # execute reconstruction
    bru.convert()
    return 0


def main(argv):
    """ request a directory containing multiple studies and the name of the destination reconstruction """
    input_dir = ''
    recon_dir = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('mass_recon -i <input_dir> -o <recon_dir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('usage: python mass_recon.py -i <input_dir> -o <recon_dir>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_dir = arg
        elif opt in ("-o", "--ofile"):
            recon_dir = arg
    
    init_dir = os.getcwd()
    cohort_study_dir = input_dir
    output_folder_name = recon_dir
    cohort_numerals = os.listdir(cohort_study_dir)
    for curr_study_idx in range(cohort_numerals.__len__()):
        curr_study_dir = os.path.join(cohort_study_dir, cohort_numerals[curr_study_idx])
        reconstruct_study(curr_study_dir, output_folder_name)
        rename_scan_dirs(curr_study_dir, output_folder_name)
    os.chdir(init_dir)
    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
