# Contributing to Simulation.jl

Thanks for contributing to Simulation.jl! Please read the below on how best to make progress on this project.

## Pull Request Process


1. Please branch from within ScottishCovidResponse/Simulation.jl 

   - Please **don't** fork this repository to your own github.com/username or other organisation.

2. Read the issues and choose/be assigned an appropriate issue. This is likely to be either a `Starter issue` or one assigned by chatting in `zulip/Simulation.jl`.

3. Name the branch `username/featurename` and have fun!

4. If referencing the issue from within a commit message then it should have the full path to the issue, `ScottishCovidResponse/SCRCIssueTracking#123`, but only in the description (not the first line of the commit) i.e.:

   - `git commit -m "pithy description" -m "ScottishCovidResponse/SCRCIssueTracking#123 and more details"` to refer to issue number 123. 
   - writing the commit message and description in an editor is also possible with `export EDITOR=vi` for example, and then just doing a `git commit`. 
   - bash also recognises line breaks in the double inverted commas

   ```bash
   $ git commit -m "first line
   >
   > second line ScottishCovidResponse/SCRCIssueTracking#123"
   ```

5. When your feature is ready to merge back into the `dev` branch create a PR against  `ScottishCovidResponse/Simulation#master` (**not** `boydorr/Simulation.jl`) and assign `ScottishCovidResponse/simulation-jl-admins` as reviewers.

6. Semver will be handled in PRs from `dev` into `master`.

7. Thanks for your time and effort!

