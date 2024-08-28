#!/usr/bin/env python

"""
Script name:  branch_PR_issue_closer.py

Goal:  To check if the newly-merged PR's commit message attempted to close an issue.
       If so, then move the associated project card to the "closed issues" column.

       Also checks if the newly-merged PR is the final PR needed to fix the issue
       for all related branches.  If so, then the issue is formally closed.

       Finally, this script also checks to see if the merged PR attempted
       to close other PRs, and does so if the merge was not to the repo's default branch.

Written by:  Jesse Nusbaumer <nusbaume@ucar.edu> - October, 2019
"""

#+++++++++++++++++++++
#Import needed modules
#+++++++++++++++++++++

import re
import sys
import argparse

from github import Github

#################
#HELPER FUNCTIONS
#################

#++++++++++++++++++++++++++++++
#Input Argument parser function
#++++++++++++++++++++++++++++++

def parse_arguments():

    """
    Parses command-line input arguments using the argparse
    python module and outputs the final argument object.
    """

    #Create parser object:
    parser = argparse.ArgumentParser(description='Close issues and pull requests specified in merged pull request.')

    #Add input arguments to be parsed:
    parser.add_argument('--access_token', metavar='<GITHUB_TOKEN>', action='store', type=str,
                        help="access token used to access GitHub API")

    parser.add_argument('--trigger_sha', metavar='<GITHUB SHA>', action='store', type=str,
                        help="Commit SHA that triggered the workflow")

    #Parse Argument inputs
    args = parser.parse_args()
    return args

#++++++++++++++++++++++++++++++++
#Script message and exit function
#++++++++++++++++++++++++++++++++

def end_script(msg):

    """
    Prints message to screen, and then exits script.
    """
    print(f"\n{msg}\n")
    print("Issue closing check has completed successfully.")
    sys.exit(0)

#############
#MAIN PROGRAM
#############

def _main_prog():

    # pylint: disable=too-many-locals
    # pylint: disable=too-many-branches
    # pylint: disable=too-many-statements

    #++++++++++++
    #Begin script
    #++++++++++++

    print("Checking if issue needs to be closed...")

    #+++++++++++++++++++++++
    #Read in input arguments
    #+++++++++++++++++++++++

    args = parse_arguments()

    #Add argument values to variables:
    token = args.access_token
    trigger_sha = args.trigger_sha

    #++++++++++++++++++++++++++++++++
    #Log-in to github API using token
    #++++++++++++++++++++++++++++++++

    ghub = Github(token)

    #+++++++++++++++++++++
    #Open ESCOMP/CAM repo
    #+++++++++++++++++++++

    cam_repo = ghub.get_repo("ESCOMP/CAM")

    #+++++++++++++++++++++++++++++
    #Get triggering commit message
    #+++++++++++++++++++++++++++++

    github_commit = cam_repo.get_commit(trigger_sha)

    commit_message = github_commit.commit.message

    #+++++++++++++++++++++++++++++++
    #Search for github PR merge text
    #+++++++++++++++++++++++++++++++

    #Compile Pull Request merge text expression:
    pr_merge_pattern = re.compile(r'Merge pull request ')

    #Search for merge text, starting at beginning of message:
    commit_msg_match = pr_merge_pattern.match(commit_message)

    #Initialize variables:
    pr_num = 0

    #Check if match exists:
    if commit_msg_match is not None:
        #If it does then pull out text immediately after message:
        post_msg_text = commit_message[commit_msg_match.end():]

        #Split text into individual words:
        post_msg_word_list = post_msg_text.split()

        #Extract first word:
        first_word = post_msg_word_list[0]

        #Print merged pr number to screen:
        print(f"Merged PR: {first_word}")

        try:
            #Try assuming the word is just a number:
            pr_num = int(first_word[1:]) #ignore "#" symbol
        except ValueError:
            #If the conversion fails, then this is likely not a real PR merge, so end the script:
            endmsg = "No Pull Request number was found in the commit message, so there is nothing for the script to do."
            end_script(endmsg)

    else:
        endmsg = "This push commit does not appear to be a merged pull request, so the script will do nothing."
        end_script(endmsg)

    #+++++++++++++++++++++++++++++++++++++
    #Check that PR has in fact been merged
    #+++++++++++++++++++++++++++++++++++++

    #Extract pull request info:
    merged_pull = cam_repo.get_pull(pr_num)

    #If pull request has not been merged, then exit script:
    if not merged_pull.merged:
        endmsg = "Pull request in commit message was not actually merged, so the script will not close anything."
        end_script(endmsg)

    #++++++++++++++++++++++++++++++++++++++++
    #Check that PR was not for default branch
    #++++++++++++++++++++++++++++++++++++++++

    #Determine default branch on repo:
    default_branch = cam_repo.default_branch

    #Extract merged branch from latest Pull request:
    merged_branch = merged_pull.base.ref

    #If PR was to default branch, then exit script (as github will handle it automatically):
    if merged_branch == default_branch:
        endmsg = "Pull request ws merged into default repo branch. Thus issue is closed automatically"
        end_script(endmsg)

    #++++++++++++++++++++++++++++++++++++++
    #Create integer list of all open issues:
    #++++++++++++++++++++++++++++++++++++++

    #Extract list of open issues from repo:
    open_repo_issues = cam_repo.get_issues(state='open')

    #Collect all open repo issues:
    open_issues = [issue.number for issue in open_repo_issues]

    #+++++++++++++++++++++++++++++++++++++++++++++
    #Create integer list of all open pull requests
    #+++++++++++++++++++++++++++++++++++++++++++++

    #Extract list of open PRs from repo:
    open_repo_pulls = cam_repo.get_pulls(state='open')

    #Collect all open pull requests:
    open_pulls = [pr.number for pr in open_repo_pulls]

    #+++++++++++++++++++++++++++++++++++++++++++++++++
    #Check if one of the keywords exists in PR message
    #+++++++++++++++++++++++++++++++++++++++++++++++++

    #Keywords are:
    #close, closes, closed
    #fix, fixes, fixed
    #resolve, resolves, resolved

    #Create regex pattern to find keywords:
    keyword_pattern = re.compile(r'(^|\s)close(\s|s\s|d\s)|(^|\s)fix(\s|es\s|ed\s)|(^|\s)resolve(\s|s\s|d\s)')

    #Extract (lower case) Pull Request message:
    pr_msg_lower = merged_pull.body.lower()

    #search for at least one keyword:
    word_matches = []
    if keyword_pattern.search(pr_msg_lower) is not None:
        #If at least one keyword is found, then determine location of every keyword instance:
        word_matches = keyword_pattern.finditer(pr_msg_lower)
    else:
        endmsg = "Pull request was merged without using any of the keywords.  Thus there are no issues to close."
        end_script(endmsg)

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Extract issue and PR numbers associated with found keywords in merged PR message
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #create issue pattern ("the number symbol {#} + a number"),
    #which ends with either a space, a comma, a period, or
    #the end of the string itself:
    issue_pattern = re.compile(r'#[0-9]+(\s|,|$)|.')

    #Create new "close" issues list:
    close_issues = []

    #Create new "closed" PR list:
    close_pulls = []

    #Search text right after keywords for possible issue numbers:
    for match in word_matches:

        #create temporary string starting at end of match:
        tmp_msg_str = pr_msg_lower[match.end():]

        #Check if first word matches issue pattern:
        if issue_pattern.match(tmp_msg_str) is not None:

            #If so, then look for an issue number immediately following
            first_word = tmp_msg_str.split()[0]

            #Extract issue number from first word:
            try:
                #First try assuming the string is just a number
                issue_num = int(first_word[1:]) #ignore "#" symbol
            except ValueError:
                #If not, then ignore last letter:
                try:
                    issue_num = int(first_word[1:-1])
                except ValueError:
                    #If ignoring the first and last letter doesn't work,
                    #then the match was likely a false positive,
                    #so set the issue number to one that will never be found:
                    issue_num = -9999

            #Check if number is actually for a PR (as opposed to an issue):
            if issue_num in open_pulls:
                #Add PR number to "close pulls" list:
                close_pulls.append(issue_num)
            elif issue_num in open_issues:
                #If in fact an issue, then add to "close issues" list:
                close_issues.append(issue_num)

    #If no issue numbers are present after any of the keywords, then exit script:
    if not close_issues and not close_pulls:
        endmsg = "No open issue or PR numbers were found in the merged PR message.  Thus there is nothing to close."
        end_script(endmsg)

    #Print list of referenced issues to screen:
    if close_issues:
        print("Issues referenced by the merged PR: "+", ".join(\
              str(issue) for issue in close_issues))

    #Print list of referenced PRs to screen:
    if close_pulls:
        print("PRs referenced by the merged PR: "+", ".join(\
              str(pull) for pull in close_pulls))

    #++++++++++++++++++++++++++++++++++++++++++++++
    #Attempt to close all referenced issues and PRs
    #++++++++++++++++++++++++++++++++++++++++++++++

    #Loop over referenced issues:
    for issue_num in close_issues:
        #Extract github issue object:
        cam_issue = cam_repo.get_issue(number=issue_num)
        #Close issue:
        cam_issue.edit(state='closed')
        print(f"Issue #{issue_num} has been closed.")

    #Loop over referenced PRs:
    for pull_num in close_pulls:
        #Extract Pull request object:
        cam_pull = cam_repo.get_pull(number=pull_num)

        #Close Pull Request:
        cam_pull.edit(state='closed')
        print(f"Pull Request #{pull_num} has been closed.")

    #++++++++++
    #End script
    #++++++++++

    print("Issue closing check has completed successfully.")

#############################################

#Run the main script program:
if __name__ == "__main__":
    _main_prog()
