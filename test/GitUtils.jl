# SPDX-License-Identifier: LGPL-3.0-or-later

module GitUtils

using Git

export is_repo_clean

"""
    is_repo_clean(repo_path; strict = false)

Determine whether the git repository at `repo_path` is clean, reporting on its state.

The short-format status (`git status -s`) of the repository is split into three
categories using the two status columns of each entry — the first column records
changes staged in the index, the second records changes in the working tree:

  - **staged**: tracked files with changes recorded in the index;
  - **unstaged**: tracked files with working-tree changes not yet staged (this
    includes files that are both staged and further modified);
  - **untracked**: files git is not tracking (shown as `??`).

By default (`strict = false`) only *unstaged* changes make the repository dirty:
a repository whose only changes are already staged and/or untracked is treated as
clean. This suits a formatting check, which should fail only when reformatting
produced changes that have not yet been staged. Pass `strict = true` to require
the repository to be *genuinely* clean, with no staged, unstaged or untracked
changes at all.

Regardless of the result, an informational message (`@info`) always reports how
many files are staged and lists any untracked files. When the repository is
considered dirty, the offending files are additionally reported via `@error`.

Return `true` if the repository is clean under the selected criterion, and `false`
otherwise.
"""
function is_repo_clean(repo_path; strict = false)
    # Short-format porcelain status: each entry is "XY path", where column X is the
    # index (staged) status and column Y is the working-tree status.
    statuses = readlines(`$(Git.git()) status -s $repo_path`)

    untracked = filter(s -> startswith(s, "??"), statuses)
    staged = filter(s -> !startswith(s, "??") && s[1] != ' ', statuses)
    unstaged = filter(s -> !startswith(s, "??") && s[2] != ' ', statuses)

    # Always report the staged count and any untracked files
    msg =
        "Repository $repo_path: $(length(staged)) file(s) staged, " *
        "$(length(untracked)) untracked"
    isempty(untracked) || (msg *= "\nUntracked files:\n" * join(untracked, "\n"))
    @info msg

    # Relaxed (default): only unstaged changes to tracked files count as dirty.
    # Strict: any staged, unstaged or untracked change counts.
    dirty = strict ? statuses : unstaged

    is_clean = isempty(dirty)
    is_clean || @error "Repository not clean:\n" * join(dirty, "\n")

    return is_clean
end

end
