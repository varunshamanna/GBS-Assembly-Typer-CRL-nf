#!/usr/bin/env bash

# Make a release.

function validate_new_version {
    local VERSION=${1}

    # get the latest tags
    git fetch --tags

    # check the version doesn't exist
    if git tag --list | grep -E -q "^${VERSION}$"; then
        echo "Version ${VERSION} already exists."
        exit 1
    fi
}

function validate_branch {
    # check branch is main
    BRANCH=$(git rev-parse --abbrev-ref HEAD)
    if [[ "$BRANCH" != "main" ]]; then
        echo 'Releases can only be made from `main` branch.'
        exit 1
    fi
}

function validate_staging_empty {
    # check no staged changes
    git diff-index --quiet HEAD --
    STAGING_STATUS=$?
    if [ "$STAGING_STATUS" != "0" ]; then
        echo "Staged changes found. Please commit them separately first."
        exit 1
    fi
}

# check args count
if [ $# -ne 1 ]; then
  echo "Usage: $0 <version (without v)>"
  exit 1
fi

VERSION=$1
validate_new_version "${VERSION}"
validate_branch
validate_staging_empty

# update version number in nextflow config file
NEXTFLOW_CONFIG="nextflow.config"
sed "s@gbs-typer-sanger-nf:.*@gbs-typer-sanger-nf:${VERSION}'@g" ${NEXTFLOW_CONFIG} > /tmp/${NEXTFLOW_CONFIG}.tmp
cat /tmp/${NEXTFLOW_CONFIG}.tmp > ${NEXTFLOW_CONFIG}
git add ${NEXTFLOW_CONFIG}
git commit -m "Updated container version number for release v${VERSION}" ${NEXTFLOW_CONFIG}

# tag and commit
echo "Creating version ${VERSION}..."
git commit -m "v${VERSION}"
git tag "v${VERSION}"

# make sure origin is up to date
echo "Pushing version ${VERSION}..."
git push origin main
git push origin "v${VERSION}"

echo "Done."
