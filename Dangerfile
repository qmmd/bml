# Check basic commit message formatting
#
# vim: syntax=ruby

# Lint the commit message.
commit_lint.check

# Ensure a clean commit history.
if git.commits.any? { |c| c.message =~ /^Merge branch/ }
    warn('Please rebase to get rid of the merge commits in this PR')
end

# Look for prose issues.
#prose.lint_files markdown_files

# Look for spelling issues.
#prose.check_spelling markdown_files

# Get information about the conflict between PRs.
conflict_checker.check_conflict

# Running plugin with reviewers count specified. Find maximum two
# reviewers.
mention.run(2, [], [])

# Prints markdown containing a random post.
the_coding_love.random

# Show the greeting.
welcome_message.greet

# Local Variables:
# mode: ruby
# End:
