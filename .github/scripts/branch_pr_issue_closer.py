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
import subprocess
import shlex
import argparse

from github import Github

#################
#HELPER FUNCTIONS
#################

#+++++++++++++++++++++++++++++++++++++++++
#Curl command needed to move project cards
#+++++++++++++++++++++++++++++++++++++++++

def  project_card_move(oa_token, column_id, card_id):

    """
    Currently pyGithub doesn't contain the methods required
    to move project cards from one column to another, so
    the unix curl command must be called directly, which is
    what this function does.

    The specific command-line call made is:

    curl -H "Authorization: token OA_token" -H \
    "Accept: application/vnd.github.inertia-preview+json" \
    -X POST -d '{"position":"top", "column_id":<column_id>}' \
    https://api.github.com/projects/columns/cards/<card_id>/moves

    """

    #create required argument strings from inputs:
    github_oa_header = ''' "Authorization: token {0}" '''.format(oa_token)
    github_url_str = '''https://api.github.com/projects/columns/cards/{0}/moves'''.format(card_id)
    json_post_inputs = ''' '{{"position":"top", "column_id":{}}}' '''.format(column_id)

    #Create curl command line string:
    curl_cmdline = '''curl -H '''+github_oa_header+''' -H "Accept: application/vnd.github.inertia-preview+json" -X POST -d '''+\
                   json_post_inputs+''' '''+github_url_str

    #Split command line string into argument list:
    curl_arg_list = shlex.split(curl_cmdline)

    #Run command using subprocess:
    subprocess.run(curl_arg_list, check=True)

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
    print("\n{}\n".format(msg))
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

    #++++++++++++++++++++
    #Open ESCOMP/CAM repo
    #++++++++++++++++++++

    #Official CAM repo:
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

    #Check if match exists:
    if commit_msg_match is not None:
        #If it does then pull out text immediately after message:
        post_msg_text = commit_message[commit_msg_match.end():]

        #Split text into individual words:
        post_msg_word_list = post_msg_text.split()

        #Extract first word:
        first_word = post_msg_word_list[0]

        #Print merged pr number to screen:
        print("Merged PR: {}".format(first_word))

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
    close_issues = list()

    #Create new "closed" PR list:
    close_pulls = list()

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

            #Check that number is actually for an issue (as opposed to a PR):
            if issue_num in open_issues:
                #Add issue number to "close issues" list:
                close_issues.append(issue_num)
            elif issue_num in open_pulls:
                #If in fact a PR, then add to PR list:
                close_pulls.append(issue_num)

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

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Determine name of project associated with merged Pull Request
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #Pull-out all projects from repo:
    projects = cam_repo.get_projects()

    #Initalize modified project name:
    proj_mod_name = None

    #Loop over all repo projects:
    for project in projects:
        #Pull-out columns from each project:
        proj_columns = project.get_columns()

        #Loop over columns:
        for column in proj_columns:

            #check if column name is "Completed Tags"
            if column.name == "Completed tags":
                #If so, then extract cards:
                cards = column.get_cards()

                #Loop over cards:
                for card in cards:
                    #Extract card content:
                    card_content = card.get_content()

                    #Next, check if card number exists and matches merged PR number:
                    if card_content is not None and card_content.number == pr_num:
                        #If so, and if Project name is None, then set string:
                        if proj_mod_name is None:
                            proj_mod_name = project.name
                            #Break out of card loop:
                            break

                        #If already set, then somehow merged PR is in two different projects,
                        #which is not what this script is expecting, so just exit:
                        endmsg = "Merged Pull Request found in two different projects, so script will do nothing."
                        end_script(endmsg)

    #Print project name associated with merged PR:
    print("merged PR project name: {}".format(proj_mod_name))

    #++++++++++++++++++++++++++++++++++++++++
    #Extract repo project "To do" card issues
    #++++++++++++++++++++++++++++++++++++++++

    #Initalize issue counting dictionary:
    proj_issues_count = dict()

    #Initalize issue id to project card id dictionary:
    proj_issue_card_ids = dict()

    #Loop over all repo projects:
    for project in projects:

        #Next, pull-out columns from each project:
        proj_columns = project.get_columns()

        #Loop over columns:
        for column in proj_columns:
            #Check if column name is "To do"
            if column.name == "To do":
                #If so, then extract cards:
                cards = column.get_cards()

               #Loop over cards:
                for card in cards:
                    #Extract card content:
                    card_content = card.get_content()

                    #Next, check if card issue number matches any of the "close" issue numbers from the PR:
                    if card_content is not None and card_content.number in close_issues:

                        #If so, then check if issue number is already in proj_issues_count:
                        if card_content.number in proj_issues_count:
                            #Add one to project issue counter:
                            proj_issues_count[card_content.number] += 1

                            #Also add issue id and card id to id dictionary used for card move, if in relevant project:
                            if project.name == proj_mod_name:
                                proj_issue_card_ids[card_content.number] = card.id

                        else:
                            #If not, then append to project issues count dictionary:
                            proj_issues_count[card_content.number] = 1

                            #Also add issue id and card id to id dictionary used for card move, if in relevant project:
                            if project.name == proj_mod_name:
                                proj_issue_card_ids[card_content.number] = card.id

            #Otherwise, check if column name matches "closed issues" column:
            elif column.name == "closed issues" and project.name == proj_mod_name:
                #If so, then save column id:
                #column_id = column.id
                column_target_id = column.id

    #If no project cards are found that match the issue, then exit script:
    if not proj_issues_count:
        endmsg = "No project cards match the issue being closed, so the script will do nothing."
        end_script(endmsg)

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Check if the number of "To-do" project cards matches the total number
    #of merged PRs for each 'close' issue.
    #
    #Then, close all issues for which project cards equals merged PRs
    #
    #If not, then simply move the project card to the relevant project's
    #"closed issues" column.
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #Loop over project issues and counts that have been "closed" by merged PR:
    for issue_num, issue_count in proj_issues_count.items():

        #If issue count is just one, then close issue:
        if issue_count == 1:
            #Extract github issue object:
            cam_issue = cam_repo.get_issue(number=issue_num)
            #Close issue:
            cam_issue.edit(state='closed')
            print("Issue #{} has been closed.".format(issue_num))
        else:
            #Extract card id from id dictionary:
            card_id = proj_issue_card_ids[issue_num]

            #Then move the card on the relevant project page to the "closed issues" column:
            project_card_move(token.strip(), column_target_id, card_id)

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Finally, close all Pull Requests in "close_pulls" list:
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++

    for pull_num in close_pulls:
        #Extract Pull request object:
        cam_pull = cam_repo.get_pull(number=pull_num)

        #Close Pull Request:
        cam_pull.edit(state='closed')
        print("Pull Request #{} has been closed.".format(pull_num))

    #++++++++++
    #End script
    #++++++++++

    print("Issue closing check has completed successfully.")

#############################################

#Run the main script program:
if __name__ == "__main__":
    _main_prog()
