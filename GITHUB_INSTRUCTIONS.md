# How to Upload This Code to GitHub

## Step-by-Step Instructions for Beginners

This guide will walk you through uploading the TCGS-SEQUENTION protocol code to GitHub. No prior experience needed!

---

## Part 1: Create a GitHub Account (Skip if you have one)

1. Go to [github.com](https://github.com)
2. Click **"Sign up"** in the top right
3. Enter your email, create a password, and choose a username
4. Verify your email address
5. Done! You now have a GitHub account

---

## Part 2: Create a New Repository

1. Log in to GitHub
2. Click the **"+"** button in the top right corner
3. Select **"New repository"**
4. Fill in:
   - **Repository name:** `sequention-protocol`
   - **Description:** `Three-Gate Falsification Protocol for Synchronous Parallel Emergence`
   - **Public** (select this option)
   - **Do NOT** check "Add a README file" (we have one already)
5. Click **"Create repository"**

You'll see a page with instructions. Keep this page open!

---

## Part 3: Download and Install Git (if not installed)

### On Windows:
1. Go to [git-scm.com/download/win](https://git-scm.com/download/win)
2. Download and run the installer
3. Click "Next" through all screens (default options are fine)
4. Restart your computer

### On Mac:
1. Open Terminal (search for "Terminal" in Spotlight)
2. Type: `git --version`
3. If prompted to install, click "Install"

### On Linux:
1. Open Terminal
2. Type: `sudo apt install git` (Ubuntu/Debian) or `sudo yum install git` (Fedora)

---

## Part 4: Upload the Files

### Method A: Using GitHub Desktop (Easiest for Beginners)

1. **Download GitHub Desktop:**
   - Go to [desktop.github.com](https://desktop.github.com)
   - Download and install it
   - Sign in with your GitHub account

2. **Clone your repository:**
   - Click "File" → "Clone repository"
   - Select your `sequention-protocol` repository
   - Choose where to save it on your computer
   - Click "Clone"

3. **Add the files:**
   - Open the folder where you cloned the repository
   - Copy ALL files from the `sequention-protocol` folder I provided
   - Paste them into the cloned folder

4. **Upload (Push) to GitHub:**
   - Go back to GitHub Desktop
   - You'll see all the new files listed
   - Type a message in the "Summary" box: `Initial upload of protocol code`
   - Click **"Commit to main"**
   - Click **"Push origin"** (or "Publish branch")

5. **Done!** Your code is now on GitHub!

---

### Method B: Using Command Line

1. **Open Terminal/Command Prompt**

2. **Navigate to the folder:**
   ```bash
   cd path/to/sequention-protocol
   ```

3. **Initialize Git:**
   ```bash
   git init
   ```

4. **Add all files:**
   ```bash
   git add .
   ```

5. **Create your first commit:**
   ```bash
   git commit -m "Initial upload of protocol code"
   ```

6. **Connect to GitHub** (copy the URL from your GitHub repository page):
   ```bash
   git remote add origin https://github.com/YOUR-USERNAME/sequention-protocol.git
   ```

7. **Upload to GitHub:**
   ```bash
   git branch -M main
   git push -u origin main
   ```

8. **Enter your credentials** when prompted

9. **Done!**

---

## Part 5: Making Future Updates

When you want to update the code:

### Using GitHub Desktop:
1. Make changes to files on your computer
2. Open GitHub Desktop
3. You'll see the changes listed
4. Type a description of what you changed
5. Click "Commit to main"
6. Click "Push origin"

### Using Command Line:
```bash
# Navigate to the folder
cd path/to/sequention-protocol

# Add changed files
git add .

# Describe your changes
git commit -m "Description of what you changed"

# Upload
git push
```

---

## Part 6: Sharing Your Repository

Your code will be available at:
```
https://github.com/YOUR-USERNAME/sequention-protocol
```

Anyone can:
- View the code
- Download it
- Clone it to their computer
- Suggest changes (via "Pull Requests")

---

## Common Problems and Solutions

### "Permission denied" error
- Make sure you're logged in
- On command line, you may need to set up SSH keys or use a personal access token

### "Repository not found" error
- Double-check the repository URL
- Make sure the repository exists on GitHub

### "Files too large" error
- GitHub has a 100MB file limit
- Remove large files or use Git LFS

### Need to undo a mistake?
- On GitHub Desktop: Click "History" tab, right-click → "Undo commit"
- On command line: `git reset HEAD~1` (undoes last commit)

---

## Quick Reference Card

| Action | GitHub Desktop | Command Line |
|--------|---------------|--------------|
| Add files | Copy files to folder | `git add .` |
| Save changes | Click "Commit" | `git commit -m "message"` |
| Upload | Click "Push" | `git push` |
| Download updates | Click "Fetch/Pull" | `git pull` |

---

## Getting Help

- **GitHub Docs:** [docs.github.com](https://docs.github.com)
- **GitHub Desktop Guide:** [docs.github.com/desktop](https://docs.github.com/en/desktop)
- **Git Cheatsheet:** [education.github.com/git-cheat-sheet-education.pdf](https://education.github.com/git-cheat-sheet-education.pdf)

---

## Final Checklist

Before sharing your repository, make sure:

- [ ] All Python files are uploaded
- [ ] README.md displays correctly on GitHub
- [ ] requirements.txt is included
- [ ] LICENSE file is included
- [ ] data/ folder with ltee_genes.json is included
- [ ] No sensitive information (passwords, API keys) is in the code

---

**You did it!** 🎉

Your TCGS-SEQUENTION protocol code is now available on GitHub for the world to see and use.
