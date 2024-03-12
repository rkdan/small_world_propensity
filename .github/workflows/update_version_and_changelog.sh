#!/bin/bash

# Increment version in pyproject.toml and __init__.py
# This is a simplistic approach; you may need a more sophisticated versioning logic.
poetry version patch

# Extract the new version number
NEW_VERSION=$(poetry version -s)

# Update __init__.py with the new version
sed -i "s/__version__ = .*/__version__ = '$NEW_VERSION'/" small_world_propensity/__init__.py

# Update CHANGELOG.md
# This is a placeholder. Consider using a tool or script to generate meaningful changelog entries.
echo -e "## $NEW_VERSION\n* Your changes here\n\n$(cat CHANGELOG.md)" > CHANGELOG.md