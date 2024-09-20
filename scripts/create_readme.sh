#!/bin/bash

common/scripts/create_task_readme --input src/api

sed -i 's#flowchart LR#flowchart TB#' README.md
