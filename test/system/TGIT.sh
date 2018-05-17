#!/bin/sh
# Test for bad git repo
# Ensures that the top-level CAM directory 
# has ".git" directory and ".gitignore" file,
# and no other git files or directories.

# Return codes in use:
# 1: Not a git repository.
# 2: Missing ".git" directory.
# 3: Missing ".gitignore" file.
# 4: More than two ".git*" files.
# 5: Error from running an external command

# Utility to check return code.
# Give it the code and an error message, and it will print stuff and exit.
check_code () {
    if [ "$1" -ne 0 ]; then
        echo "Error: return code from command was $1"
        echo "$2"
        exit 5
    fi
}

# Little utility for finding absolute path to a directory.
get_dir_abspath () {
    echo $(cd $1 && pwd)
}

# Set CAM top-level directory:
cam_top_dir=$(get_dir_abspath ${CAM_SCRIPTDIR}/../..)

# Initialize error variable:
rc=0

# Check to make sure that the top level directory is a git repo:
top_dir_gitinfo=$(git -C "${cam_top_dir}" rev-parse)

# Save error code:
gitinfo_error=$?

# Print error if directory is not a git repository:
if [ "${gitinfo_error}" -ne 0 ]; then
    cat <<EOF
The top directory does not appear to be a git repository. Make sure:
1) Git is present on your system and PATH (i.e. "git" exists as a command in your shell).
2) You are running these tests on the correct directory, and
3) The sub-directory ".git" exists in the directory (if not then copy it from a proper CAM git repo).
EOF
    rc=1
else # If directory is a git repoistory, then run other tests:

    # Check for ".git" directory (i.e., make sure git repo isn't "bare"):
    if [ ! -d "${cam_top_dir}/.git" ]; then
        cat <<EOF
The ".git" directory is missing from the CAM git repo.  Was this repo cloned, copied, or
modified incorrectly?  If so then copy the .git directory from a standard CAM git repo.
EOF
        rc=2
    fi

    # Check for missing ".gitignore" file.
    if [ ! -f "${cam_top_dir}/.gitignore" ]; then
        cat <<EOF
The ".gitignore" file is missing from the CAM git repo.  Was this repo cloned, copied, or
modified incorrectly?  If so then copy the .gitignore file from a standard CAM git repo.
EOF
        rc=3
    fi   

    # Check if there are more ".git*" files or directories than just ".git" or ".gitignore".
    git_file_num=$(find "${cam_top_dir}" -maxdepth 1 -name '.git*' | wc -l)

    check_code "$?" "Problem running 'find' command for multi-git file check."

    if [ "${git_file_num}" -gt 2 ]; then
        cat <<EOF
More than two ".git*" files or sub-directories present in this CAM git repo.
EOF
        rc=4
    fi

fi # git-repo check

# If all tests pass, let user know:
if [ "$rc" -eq 0 ]; then
    echo "No problems found in the git repostiory ($cam_top_dir)."
fi

# Exit script with final error value:
exit $rc
