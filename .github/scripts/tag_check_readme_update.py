#!/usr/bin/env python

"""
Script name:  tag_check_readme_update.py

Goal:  To determine if the pushed tag was configured properly.
       If not then the action will "fail", notifying the admins.

       Also, if the pushed tag was configured correctly and was
       for the "development" branch, then this script will update
       the README.md file on the default branch to display the new
       development tag.

Written by:  Jesse Nusbaumer <nusbaume@ucar.edu> - February, 2020
"""

#+++++++++++++++++++++
#Import needed modules
#+++++++++++++++++++++

import re
import sys
import argparse
import base64
import smtplib
import ssl

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

    parser.add_argument('--cam_email_pword', metavar='<CAM PASSWORD>', action='store', type=str,
                        help="Password for cam.engineering.works email")

    #Parse Argument inputs
    args = parser.parse_args()
    return args

#++++++++++++++++++++++
#Email message function
#++++++++++++++++++++++

def email_msg(receiver_list, pword, title, msg):

    """
    Send email message to receiver list via
    the "cam.engineering.works@gmail.com" email.
    """

    #Combine email title and message,
    #so that email is properly formatted:
    msg_email = "Subject: "+title+"\n\n"+msg

    #create email encryption:
    #Note that gmail automatically converts SSL to TLS
    context = ssl.create_default_context()

    #Open email server (gmail always uses port 465):
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
        #Login to cam.engineering.works:
        server.login("cam.engineering.works@gmail.com", pword)
        #Loop over receiver emails:
        for receiver in receiver_list:
            #Send email
            server.sendmail("cam.engineering.works@gmail.com", receiver,
                            msg_email)

#++++++++++++++++++++++++++++++++
#Script message and exit function
#++++++++++++++++++++++++++++++++

def end_script(msg):

    """
    Prints message to screen, and then exits script.
    """
    print("\n{}\n".format(msg))
    print("README tag update script has completed successfully.")
    sys.exit(0)

#+++++++++++++++++++++++++++++++++++++++++++
#Script message and exit function with error
#+++++++++++++++++++++++++++++++++++++++++++

def end_script_fail(cam_pword, msg):

    """
    Prints message to screen, and then raises an exception
    in order to generate a failure notification in Github.

    This function will also send an email to the CAM SEs
    notifying them that something has gone wrong.
    """

    receiver_list = ["nusbaume@ucar.edu", "cacraig@ucar.edu", "goldy@ucar.edu", "courtneyp@ucar.edu"]
    email_msg(receiver_list, cam_pword, "Tag README update failure!", msg)

    print("\n{}\n".format(msg))
    print("README tag update script was un-successful.")
    sys.exit(1)


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

    print("Checking if README file needs to be updated...")

    #+++++++++++++++++++++++
    #Read in input arguments
    #+++++++++++++++++++++++

    args = parse_arguments()

    #Add argument values to variables:
    token = args.access_token
    trigger_sha = args.trigger_sha
    cam_pword = args.cam_email_pword

    #++++++++++++++++++++++++++++++++
    #Log-in to github API using token
    #++++++++++++++++++++++++++++++++

    ghub = Github(token)

    #++++++++++++++++++++
    #Open ESCOMP/CAM repo
    #++++++++++++++++++++

    #Official CAM repo:
    cam_repo = ghub.get_repo("ESCOMP/CAM")

    #+++++++++++++++++++++
    #Get list of repo tags
    #+++++++++++++++++++++

    cam_tags = cam_repo.get_tags()

    #+++++++++++++++++++++++++++++++++++++++++++
    #Search for tag with sha that matches commit
    #+++++++++++++++++++++++++++++++++++++++++++

    tag_name = None

    for cam_tag in cam_tags:
        if cam_tag.commit.sha == trigger_sha:
            #If matching tag is found, then extract
            #tag name and tag commit message:
            tag_name = cam_tag.name
            tag_commit = cam_repo.get_commit(trigger_sha)

            #End tag loop:
            break

    #+++++++++++++++++++++++++++++++++++
    #If no tag matches, then exit script
    #+++++++++++++++++++++++++++++++++++

    if not tag_name:
        endmsg = "No tag was created by this push, so there is nothing to do."
        end_script(endmsg)
    else:
        print("Script found tag name of '{}'".format(tag_name))

    #+++++++++++++++++++++++++++++++
    #Search for associated PR number
    #+++++++++++++++++++++++++++++++

    #Extract tag commit message:
    commit_message = tag_commit.commit.message

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
        print("Merged PR associated with tag: {}".format(first_word))

        try:
            #Try assuming the word is just a number:
            pr_num = int(first_word[1:]) #ignore "#" symbol
        except ValueError:
            #If the conversion fails, then this is likely not a real PR merge, so end the script:
            endmsg = "No Pull Request number was found in the tagged commit message.\n"
            endmsg += "This is a non-normal tag commit, and so would be worth double-checking."
            end_script_fail(cam_pword, endmsg)

    else:
        endmsg = "No Pull Request merges were found in the tagged commit message.\n"
        endmsg += "This is likely a non-normal tag commit, and so would be worth double-checking."
        end_script_fail(cam_pword, endmsg)

    #+++++++++++++++++++++++++++++++++++++
    #Check if tagged PR was in fact merged
    #+++++++++++++++++++++++++++++++++++++

    #Extract pull request info:
    merged_pull = cam_repo.get_pull(pr_num)

    #If pull request has not been merged, then exit script:
    if not merged_pull.merged:
        endmsg = "Pull request in tagged commit message was not actually merged.\n"
        endmsg += "This is likely a non-normal tag commit, and so would be worth double-checking."
        end_script_fail(cam_pword, endmsg)

    #++++++++++++++++++++++++++++++++++++++++++
    #Check if PR was for the development branch
    #++++++++++++++++++++++++++++++++++++++++++

    #Extract merged branch from latest Pull request:
    merged_branch = merged_pull.base.ref

    print("Script found branch name of '{}'".format(merged_branch))

    #If PR is not to the development branch, then exit script:
    if merged_branch != "cam_development":
        endmsg = "Tagged PR merged into non-development branch. No further action will thus be taken."
        end_script(endmsg)

    #++++++++++++++++++++++
    #Extract ChangeLog file
    #++++++++++++++++++++++

    #Extract contents of "doc" directory:
    doc_files = cam_repo.get_dir_contents("doc", ref=merged_branch)

    #Extract SHA for ChangeLog:
    changelog_sha = None
    for doc_file in doc_files:
        if doc_file.name == "ChangeLog":
            changelog_sha = doc_file.sha
            break

    #If sha can't be found, then it means the ChangeLog is missing,
    #so throw an error:
    if not changelog_sha:
        endmsg = "ChangeLog SHA can't be found, which means no doc/ChangeLog file exists.\n"
        endmsg += "Please create a new PR ASAP that adds the ChangeLog back to the development branch."
        end_script_fail(cam_pword, endmsg)

    #The ChangeLog is too large to download using the standard API,
    #so we need to use the Github data API:
    changelog_git_blob = cam_repo.get_git_blob(changelog_sha)

    #Convert changelog blob to unicode string:
    changelog = base64.b64decode(changelog_git_blob.content).decode("UTF-8")

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #Search for tag in ChangeLog, to ensure log was properly updated
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    #Make ChangeLog all lowercase, so that searching is case independent:
    changelog = changelog.lower()

    #Find first occurence of "Tag name" string:
    tag_log_idx = changelog.find("tag name")

    if tag_log_idx == -1:
        endmsg = "No 'tag name' entry found in the ChangeLog!\n"
        endmsg += "The ChangeLog file has likely been corrupted."
        end_script_fail(cam_pword, endmsg)
    else:
        #Find first occurence of "Originator", to set tag search end:
        orig_log_idx = changelog.find("originator")
        if orig_log_idx == -1:
            endmsg = "No 'Originator' entry found in the ChangeLog!\n"
            endmsg += "The ChangeLog file has likely been corrupted."
            end_script_fail(cam_pword, endmsg)

    #Search for tag string inbetween "tag name" and "originator":
    tag_found = changelog[tag_log_idx:orig_log_idx].find(tag_name)

    if tag_found == -1:
        endmsg = "Newly-created tag '{}' not the first tag listed at top of ChangeLog!\n"
        endmsg += "Please update the ChangeLog so that the newest tag is listed first."
        end_script_fail(cam_pword, endmsg.format(tag_name))

    #+++++++++++++++++++++++++++++++++
    #Extract README file contents/text
    #+++++++++++++++++++++++++++++++++

    #Grab README file object:
    #Note: This command always uses the default branch unless
    #      specifically told to do otherwise.
    readme_obj = cam_repo.get_contents("README.md")

    #Extract README content (as a bytestring):
    readme_content = readme_obj.decoded_content

    #Convert bytestring to unicode (regular) string:
    readme_text = readme_content.decode('UTF-8')

    #+++++++++++++++++++++++++++++++++
    #Upate README file development tag
    #+++++++++++++++++++++++++++++++++

    #Compile development tag text expression:
    dev_tag_pattern = re.compile(r'cam\d+_\d+_\d+')

    #Search for tag text in README file:
    dev_tag_match = dev_tag_pattern.search(readme_text)

    #End script if nothing is found:
    if not dev_tag_match:
        endmsg = "No text matches expected development tag pattern in README.md file, "
        endmsg += "so no further action will be taken."
        end_script(endmsg)

    #Extract start and end indices of matched tag text:
    dev_tag_idx = dev_tag_match.span()

    #Add new tag to README text:
    new_content = readme_text[:dev_tag_idx[0]] + tag_name + readme_text[dev_tag_idx[1]:]

    #+++++++++++++++++++++++++++++++++++++++++++++++
    #Push updated README file back to default branch
    #+++++++++++++++++++++++++++++++++++++++++++++++

    cam_repo.update_file("README.md", "Auto-update development tag.", new_content, readme_obj.sha)

    #++++++++++
    #End script
    #++++++++++

    endmsg = "cam_development tag has been successfully updated in README file."

    #Google will turn off "non-approved" SMTP access for a gmail account
    #if not used regularly.  So send an email every time this action is
    #triggered and updates the README file to ensure that google keeps the access open:
    email_msg(["nusbaume@ucar.edu"], cam_pword, "Tag README update success!", endmsg)

    #Print end message to Github Action stdout:
    print(endmsg)

#############################################

#Run the main script program:
if __name__ == "__main__":
    _main_prog()
