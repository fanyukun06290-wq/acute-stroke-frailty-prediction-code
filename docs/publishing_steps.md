# Publishing steps

## Recommended release path

Use GitHub for the readable working repository and Zenodo for the citable archived version with a DOI. Publish only the cleaned code tree; do not upload `data/` or `results/` contents.

## Before the first commit

1. The release uses the MIT license for the analysis code only. The corresponding author should still confirm that the institutional policy permits this license; the license does not cover MIMIC-IV or eICU-CRD data.
2. Replace `[INSERT PUBLIC REPOSITORY URL]` in `CITATION.cff` and `docs/code_availability_statement.md` only after the repository URL is known.
3. Check the staged file list and confirm that no CSV, RDS, database export, password, Navicat screenshot, or local absolute path is present.

## GitHub using GitHub Desktop

1. Create an empty GitHub repository. Do not initialize it with a README, license, or `.gitignore`, because those files are already included here.
2. In GitHub Desktop, choose **File > New repository**, set the local path to this release folder, and initialize the repository there. If the folder has already been initialized with Git, use **Add > Add existing repository** instead.
3. Review the changed files. Commit the cleaned tree with a message such as `Initial reproducible analysis code release`.
4. Choose **Publish repository**, select the intended GitHub account or organization, and set the repository visibility to **Public** only after the final privacy check.

## GitHub using the command line

From this release folder, after creating the empty GitHub repository, replace the placeholder URL and run:

```text
git init -b main
git add README.md .gitignore CITATION.cff LICENSE sql R notebooks docs data/README.md results/README.md
git diff --cached --name-only
git commit -m "Initial reproducible analysis code release"
git remote add origin https://github.com/<USERNAME>/acute-stroke-frailty-prediction-code.git
git push -u origin main
```

The explicit `git add` list is intentional: it prevents accidentally staging local data or result files. GitHub's official guide also warns not to push passwords or API keys.

## Zenodo DOI

1. Sign in to Zenodo and connect the GitHub account.
2. Enable this GitHub repository in the Zenodo GitHub integration.
3. Create a GitHub release, preferably `v1.0.0`, after the final code review.
4. Wait for Zenodo to archive the release and copy the DOI record URL.
5. Put the DOI and repository URL into `CITATION.cff`, `docs/code_availability_statement.md`, and the manuscript's Code Availability section. Do not change the archived files after publication; create a new version if the code changes.

Official instructions:

- GitHub: https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github
- Zenodo repository integration: https://help.zenodo.org/docs/github/enable-repository/
- Zenodo GitHub release archiving: https://help.zenodo.org/docs/github/archive-software/github-upload/
