#!/usr/bin/env python

"""
Script name:  branch_PR_issue_closer.py

Goal:  To check if the newly-merged PR's commit message attempted to close an issue.
       If so, then move the associated project card to the "closed issues" column.

       Also checks if the newly-merged PR is the final PR needed to fix the issue
       for all related branches.  If so, then the issue is formally closed.

       Finally, this script also checks to see if the merged PR attempted
       to close other PRs, and does so if the merge was not to the repo's default branch.
       Any PRs referenced this way (including PRs that are already closed or merged)
       are also searched (recursively) for additional issues and PRs to close.

Written by:  Jesse Nusbaumer <nusbaume@ucar.edu> - October, 2019
"""

#+++++++++++++++++++++
#Import needed modules
#+++++++++++++++++++++

import re
import sys
import argparse

from github import Github, GithubException

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

#++++++++++++++++++++++++++++++++++++
#Function to add a list of referenced Github
#numbers to the official PR search list
#++++++++++++++++++++++++++++++++++++

def add_closed_pulls(nums, cam_repo, checked_nums, searched_pulls, search_pulls):

    """
    Add any numbers that are closed PRs to the search list.

    nums           -> list of Github issue/PR numbers to check
    cam_repo       -> PyGithub repository object
    checked_nums   -> set of numbers already checked (modified in place)
    searched_pulls -> set of PR numbers already searched
    search_pulls   -> list of PR numbers to search (modified in place)
    """

    for num in nums:

        #Ignore issues/PRs that have already been checked:
        if num in checked_nums or num in searched_pulls:
            continue

        #Add number to set of checked Github issue/PR numbers:
        checked_nums.add(num)

        try:
            #The issue object is only a PR if it has "pull_request" data:
            if cam_repo.get_issue(number=num).pull_request is not None:
                search_pulls.append(num)
        except GithubException:
            #Number doesn't exist in repo (e.g. a typo), so ignore it:
            print(f"Referenced number #{num} was not found in the repo, so it will be ignored.")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Function to find issue and PR numbers marked for closing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def find_close_numbers(pr_body, open_issues, open_pulls):

    """
    Searches the text of a pull request body for any of the GitHub
    closing keywords, and returns three lists of the numbers that
    immediately follow those keywords: the open issue numbers, the
    open PR numbers, and all other (i.e. closed issue or PR) numbers.
    """

    #Initialize output lists:
    close_issues = []
    close_pulls = []
    other_nums = []

    #Keywords are:
    #close, closes, closed
    #fix, fixes, fixed
    #resolve, resolves, resolved

    #Create regex pattern to find keywords:
    keyword_pattern = re.compile(r'(^|\s)close(\s|s\s|d\s)|(^|\s)fix(\s|es\s|ed\s)|(^|\s)resolve(\s|s\s|d\s)')

    #Extract (lower case) Pull Request message
    #(the body is "None" if the PR has no description):
    pr_msg_lower = (pr_body or "").lower()

    #create issue pattern ("the number symbol {#} + a number"),
    #which ends with either a space, a comma, a period, or
    #the end of the string itself:
    issue_pattern = re.compile(r'#[0-9]+(\s|,|$)|.')

    #Search text right after keywords for possible issue numbers:
    for match in keyword_pattern.finditer(pr_msg_lower):

        #create temporary string starting at end of match:
        tmp_msg_str = pr_msg_lower[match.end():]

        #Check if first word matches issue pattern:
        if issue_pattern.match(tmp_msg_str) is not None:

            #If so, then look for an issue number immediately following,
            #skipping when there are no words to extract a number from:
            tmp_msg_words = tmp_msg_str.split()
            if not tmp_msg_words:
                continue
            first_word = tmp_msg_words[0]

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
            elif issue_num > 0:
                #If neither, then it may be an already-closed PR, so
                #add to "other numbers" list:
                other_nums.append(issue_num)

    return close_issues, close_pulls, other_nums

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

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Extract issue and PR numbers associated with keywords in the merged PR message
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    close_issues, close_pulls, other_nums = find_close_numbers(merged_pull.body,
                                                               open_issues, open_pulls)

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Recursively search all referenced PRs (open, closed, or merged) for
    #additional issues and PRs to close
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #Keep track of PRs that have already been searched,
    #in order to prevent infinite loops from circular references:
    searched_pulls = {pr_num}

    #Keep track of numbers that have already been checked
    #for whether they are closed PRs, to avoid repeated API calls:
    checked_nums = set()

    #Create list of PRs to search, starting with the open PRs being closed.
    #Note that this list grows as new PRs are found, so newly-added PRs
    #will also be searched by the loop below:
    search_pulls = list(close_pulls)

    #Check if non-open Github numbers are actually PRs:
    add_closed_pulls(other_nums, cam_repo, checked_nums, searched_pulls, search_pulls)

    for pull_num in search_pulls:

        #Skip PRs that have already been searched:
        if pull_num in searched_pulls:
            continue
        searched_pulls.add(pull_num)

        #Extract Pull request object:
        ref_pull = cam_repo.get_pull(number=pull_num)

        #Search referenced PR message for keywords:
        ref_issues, ref_pulls, ref_others = find_close_numbers(ref_pull.body,
                                                               open_issues, open_pulls)

        #Add any newly-found issues and open PRs to the "close" lists:
        close_issues.extend(issue for issue in ref_issues if issue not in close_issues)
        close_pulls.extend(pull for pull in ref_pulls if pull not in close_pulls)

        #Add newly-found open and closed PRs to the search list:
        search_pulls.extend(ref_pulls)
        add_closed_pulls(ref_others, cam_repo, checked_nums, searched_pulls, search_pulls)

    #Remove the merged PR itself, in case it was referenced by one of the other PRs:
    close_pulls = [pull for pull in close_pulls if pull != pr_num]

    #Remove any duplicate issue and PR numbers, while preserving order:
    close_issues = list(dict.fromkeys(close_issues))
    close_pulls = list(dict.fromkeys(close_pulls))

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
