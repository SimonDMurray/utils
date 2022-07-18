#!/bin/bash


JIRA_USERNAME=$1
JIRA_PASSWORD=$2
JIRA_URL=$3

for toplevel in TIC-*; do
  cd $toplevel
  wd=`pwd`
  for sublevel in tic-*; do 
    ACTIONS=$wd/$sublevel/actions
    JIRA_TICKET_KEY="${sublevel^^}"
    echo $JIRA_TICKET_KEY
    TICKET_JSON=$(curl -s -u "${JIRA_USERNAME}:${JIRA_PASSWORD}" \
      -H "X-Atlassian-Token: nocheck" \
      ${JIRA_URL}/rest/api/2/issue/${JIRA_TICKET_KEY})
    ASSIGNEE=$( jq -r '.fields.assignee.name' <<< $TICKET_JSON ) 
    echo "Adding file ${ACTIONS}/ASSIGNEE to ${sublevel}"
    echo $ASSIGNEE > $ACTIONS/ASSIGNEE;
  done;
  cd ../
done
