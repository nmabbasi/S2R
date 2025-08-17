# Shell2R Deployment Guide

This guide will help you deploy your Shell2R website to GitHub Pages with a custom domain.

## 📋 Prerequisites

- GitHub account
- Your domain name (e.g., shell2r.com)
- Access to your domain's DNS settings

## 🚀 Step-by-Step Deployment

### Step 1: Prepare Your Repository

1. **Create a new GitHub repository**:
   - Go to [GitHub](https://github.com) and sign in
   - Click "New repository"
   - Name it something like `shell2r` or `Shell2SingleR`
   - Make it **Public** (required for free GitHub Pages)
   - Don't initialize with README, .gitignore, or license

2. **Upload your website files**:
   - Extract the `shell2r-simple-website.tar.gz` file
   - Upload ALL contents of the `shell2r-simple` folder to your repository root
   - Your repository should contain:
     ```
     index.html
     style.css
     script.js
     pages/
     lessons/
     images/
     README.md
     DEPLOYMENT_GUIDE.md
     TUTORIAL_GUIDE.md
     ```

### Step 2: Enable GitHub Pages

1. **Go to your repository settings**:
   - Click the "Settings" tab in your repository
   - Scroll down to "Pages" in the left sidebar

2. **Configure GitHub Pages**:
   - **Source**: Select "Deploy from a branch"
   - **Branch**: Select "main" (or "master" if that's your default)
   - **Folder**: Select "/ (root)"
   - Click "Save"

3. **Wait for deployment**:
   - GitHub will start building your site
   - This usually takes 1-5 minutes
   - You'll see a green checkmark when it's ready

### Step 3: Configure Custom Domain (shell2r.com)

1. **Add custom domain in GitHub**:
   - Still in Pages settings, find "Custom domain"
   - Enter: `shell2r.com`
   - Click "Save"
   - GitHub will create a `CNAME` file in your repository

2. **Configure DNS settings**:
   - Go to your domain registrar (where you bought shell2r.com)
   - Find DNS settings or DNS management
   - Add these A records:

   ```
   Type: A
   Name: @ (or leave blank for root domain)
   Value: 185.199.108.153

   Type: A
   Name: @ (or leave blank for root domain)
   Value: 185.199.109.153

   Type: A
   Name: @ (or leave blank for root domain)
   Value: 185.199.110.153

   Type: A
   Name: @ (or leave blank for root domain)
   Value: 185.199.111.153
   ```

3. **Wait for DNS propagation**:
   - DNS changes can take 24-48 hours to fully propagate
   - You can check status using online DNS checkers
   - Your site should be accessible at `https://shell2r.com` once ready

### Step 4: Enable HTTPS

1. **Enable HTTPS enforcement**:
   - Back in GitHub Pages settings
   - Check "Enforce HTTPS" (this option appears after DNS is configured)
   - This ensures your site always uses secure connections

## 🔧 Troubleshooting

### Common Issues

**"The custom domain shell2r.com is not properly formatted"**
- Make sure the `CNAME` file contains only `shell2r.com` (no http:// or www.)
- Wait a few minutes and try again
- Clear your browser cache

**Tutorials not loading on GitHub Pages**
- Check that all files are in the repository root
- Verify tutorial files are in the `lessons/` folder
- Check browser console for error messages
- Ensure file names match exactly in `script.js`

**Site not updating after changes**
- GitHub Pages can take 5-10 minutes to update
- Check the Actions tab for build status
- Hard refresh your browser (Ctrl+F5 or Cmd+Shift+R)

**DNS not resolving**
- DNS changes take time (up to 48 hours)
- Use `nslookup shell2r.com` to check DNS status
- Contact your domain registrar if issues persist

### Verification Steps

1. **Check GitHub Pages build**:
   - Go to your repository's "Actions" tab
   - Look for "pages build and deployment" workflows
   - Green checkmark = successful, red X = failed

2. **Test your website**:
   - Visit your GitHub Pages URL (usually `username.github.io/repository-name`)
   - Test all navigation links
   - Verify tutorials load correctly
   - Check mobile responsiveness

3. **Verify custom domain**:
   - Visit `https://shell2r.com`
   - Check that HTTPS is working
   - Test all pages and functionality

## 📝 Updating Your Website

### Adding New Tutorials

1. **Create new tutorial file** in `lessons/` folder
2. **Update `script.js`** to include the new file
3. **Commit and push** changes to GitHub
4. **Wait 5-10 minutes** for GitHub Pages to update

### Modifying Content

1. **Edit files** directly on GitHub or locally
2. **Commit changes** to your repository
3. **Push to main branch**
4. **Changes appear automatically** on your live site

### Updating Static Pages

1. **Edit files** in the `pages/` folder
2. **Maintain consistent styling** with existing pages
3. **Test locally** before pushing to GitHub
4. **Commit and push** changes

## 🎯 Best Practices

### Performance
- Optimize images before uploading
- Keep tutorial files reasonably sized
- Test loading speed regularly

### SEO
- Update meta descriptions in HTML files
- Use descriptive file names
- Add alt text to images

### Maintenance
- Regularly check for broken links
- Update tutorial content as needed
- Monitor Google Search Console for issues

### Security
- Keep legal pages updated
- Monitor for any security advisories
- Use HTTPS everywhere

## 📞 Getting Help

If you encounter issues:

1. **Check GitHub Pages documentation**: [docs.github.com/pages](https://docs.github.com/pages)
2. **Review GitHub Actions logs** for build errors
3. **Test locally first** to isolate issues
4. **Check browser console** for JavaScript errors
5. **Verify DNS settings** with your domain registrar

## 🎉 Success!

Once everything is set up:
- Your website will be live at `https://shell2r.com`
- New tutorials automatically appear when you add them
- The site is optimized for Google AdSense
- All legal pages are included for compliance
- Mobile-responsive design works on all devices

Your Shell2R educational platform is now ready to help others learn bioinformatics!

