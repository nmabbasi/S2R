# Shell2R Deployment Guide

Complete guide for deploying your Shell2R website to GitHub Pages with custom domain support.

## 📋 Prerequisites

Before starting deployment:

- [x] Website files are ready and tested locally
- [x] GitHub account created
- [x] Custom domain purchased (optional)
- [x] All content reviewed and finalized

## 🚀 Step-by-Step Deployment

### Step 1: Create GitHub Repository

1. **Go to GitHub.com** and sign in to your account

2. **Create a new repository**:
   - Click the "+" icon in the top right
   - Select "New repository"
   - Repository name: `shell2r` (or your preferred name)
   - Description: "Bioinformatics learning platform"
   - Make it **Public** (required for free GitHub Pages)
   - Don't initialize with README (we have our own files)
   - Click "Create repository"

### Step 2: Upload Website Files

#### Option A: Using Git Command Line

```bash
# Navigate to your website directory
cd shell2r-simple

# Initialize git repository
git init

# Add remote origin (replace 'yourusername' with your GitHub username)
git remote add origin https://github.com/yourusername/shell2r.git

# Add all files
git add .

# Commit files
git commit -m "Initial commit: Shell2R bioinformatics website"

# Push to GitHub
git push -u origin main
```

#### Option B: Using GitHub Web Interface

1. **Drag and drop files**:
   - Go to your empty repository on GitHub
   - Click "uploading an existing file"
   - Drag all files from `shell2r-simple/` folder
   - Write commit message: "Initial commit: Shell2R website"
   - Click "Commit changes"

#### Option C: Using GitHub Desktop

1. **Clone repository**:
   - Open GitHub Desktop
   - Clone your repository locally
   - Copy all files from `shell2r-simple/` to the cloned folder
   - Commit and push changes

### Step 3: Enable GitHub Pages

1. **Go to repository settings**:
   - Navigate to your repository on GitHub
   - Click the "Settings" tab (near the top right)

2. **Configure Pages**:
   - Scroll down to "Pages" in the left sidebar
   - Under "Source", select "Deploy from a branch"
   - Branch: Select "main" (or "master" if that's your default)
   - Folder: Select "/ (root)"
   - Click "Save"

3. **Wait for deployment**:
   - GitHub will show a message: "Your site is ready to be published"
   - Initial deployment takes 5-10 minutes
   - You'll get a URL like: `https://yourusername.github.io/shell2r/`

### Step 4: Test Your Deployed Site

1. **Visit your site**:
   - Click the provided URL or visit `https://yourusername.github.io/shell2r/`
   - Verify all pages load correctly
   - Test navigation, search, and tutorial loading
   - Check mobile responsiveness

2. **Common issues to check**:
   - All images display correctly
   - Tutorials load and display properly
   - Search functionality works
   - All navigation links work
   - AdSense placeholders are visible

## 🌐 Custom Domain Setup (Optional)

If you want to use your own domain (e.g., `shell2r.com`):

### Step 1: Add CNAME File

Create a file named `CNAME` (no extension) in your repository root:

```bash
# Content of CNAME file (replace with your domain)
shell2r.com
```

Commit and push this file to GitHub.

### Step 2: Configure DNS

Log in to your domain registrar (GoDaddy, Namecheap, etc.) and add these DNS records:

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

**For www subdomain** (optional):
```
Type: CNAME
Name: www
Value: yourusername.github.io
```

### Step 3: Configure Custom Domain in GitHub

1. **Go to repository Settings > Pages**
2. **Under "Custom domain"**:
   - Enter your domain: `shell2r.com`
   - Click "Save"
3. **Wait for DNS propagation** (24-48 hours)
4. **Enable HTTPS**:
   - Once DNS propagates, check "Enforce HTTPS"
   - This ensures secure connections to your site

### Step 4: Verify Custom Domain

1. **Check DNS propagation**:
   - Use tools like `whatsmydns.net` to verify DNS changes
   - All regions should show your GitHub Pages IP addresses

2. **Test your domain**:
   - Visit `https://yourdomain.com`
   - Verify SSL certificate is working (green lock icon)
   - Test all functionality

## 🔧 Advanced Configuration

### Environment-Specific Settings

You might want different settings for local development vs. production:

```javascript
// In script.js, detect environment
const isProduction = window.location.hostname !== 'localhost' && 
                    window.location.hostname !== '127.0.0.1';

// Use different API endpoints or configurations
const config = {
    baseUrl: isProduction ? 'https://shell2r.com' : 'http://localhost:8000',
    // Add other environment-specific settings
};
```

### Google Analytics Integration

Add Google Analytics to track visitors:

```html
<!-- Add to <head> section of all HTML files -->
<!-- Google tag (gtag.js) -->
<script async src="https://www.googletagmanager.com/gtag/js?id=GA_MEASUREMENT_ID"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'GA_MEASUREMENT_ID');
</script>
```

### Search Console Setup

1. **Add property to Google Search Console**:
   - Go to `search.google.com/search-console`
   - Add your domain as a property
   - Verify ownership using HTML file method

2. **Submit sitemap**:
   - Create a simple sitemap.xml file
   - Submit it to Search Console for better indexing

## 🔄 Updating Your Website

### Adding New Tutorials

1. **Create new markdown file** in `lessons/` directory
2. **Add images** to `images/` directory if needed
3. **Commit and push changes**:
   ```bash
   git add .
   git commit -m "Add new tutorial: [Tutorial Name]"
   git push
   ```
4. **GitHub Pages automatically rebuilds** (takes 2-5 minutes)

### Updating Content

1. **Edit existing files** locally
2. **Test changes** using local server
3. **Commit and push**:
   ```bash
   git add .
   git commit -m "Update: [description of changes]"
   git push
   ```

### Managing Versions

Use git tags for important releases:

```bash
# Tag a release
git tag -a v1.0 -m "Initial release"
git push origin v1.0

# List all tags
git tag -l
```

## 🛡️ Security Best Practices

### Repository Security

1. **Never commit sensitive data**:
   - API keys
   - Passwords
   - Personal information

2. **Use .gitignore** for local files:
   ```
   # .gitignore
   .DS_Store
   Thumbs.db
   *.log
   node_modules/
   .env
   ```

3. **Enable branch protection** (optional):
   - Go to Settings > Branches
   - Add rule for main branch
   - Require pull request reviews

### Website Security

1. **HTTPS enforcement**:
   - Always enable "Enforce HTTPS" in GitHub Pages settings
   - Update all internal links to use HTTPS

2. **Content Security Policy**:
   - Consider adding CSP headers for additional security
   - Be careful with external CDN resources

## 📊 Monitoring and Analytics

### GitHub Pages Analytics

Monitor your deployment:

1. **Repository Insights**:
   - View traffic statistics
   - Monitor clone/download activity
   - Track popular content

2. **Actions tab**:
   - Monitor deployment status
   - View build logs if issues occur

### Website Performance

1. **Google PageSpeed Insights**:
   - Test your site regularly
   - Optimize based on recommendations

2. **Lighthouse audits**:
   - Run audits for performance, accessibility, SEO
   - Address any issues found

## 🐛 Troubleshooting Deployment

### Common Issues and Solutions

#### 1. Site Not Loading After Deployment

**Problem**: 404 error or blank page

**Solutions**:
- Verify `index.html` is in repository root
- Check GitHub Pages settings are correct
- Wait 10-15 minutes for initial deployment
- Check repository is public

#### 2. Custom Domain Not Working

**Problem**: Domain shows 404 or doesn't resolve

**Solutions**:
- Verify CNAME file contains only your domain name
- Check DNS records are correct
- Wait for DNS propagation (up to 48 hours)
- Use DNS checker tools to verify propagation

#### 3. HTTPS Not Available

**Problem**: "Enforce HTTPS" option is grayed out

**Solutions**:
- Wait for DNS propagation to complete
- Verify custom domain is properly configured
- Check that CNAME record points to GitHub Pages

#### 4. Tutorials Not Loading

**Problem**: "Loading tutorials..." persists on live site

**Solutions**:
- Verify markdown files are in `lessons/` directory
- Check file names match JavaScript expectations
- Ensure files are committed and pushed to repository
- Check browser console for errors

#### 5. Images Not Displaying

**Problem**: Broken image links on live site

**Solutions**:
- Verify images are in `images/` directory
- Check image file names and paths in markdown
- Ensure images are committed to repository
- Use relative paths (not absolute)

### Getting Help

If you encounter issues:

1. **Check GitHub Status**: `githubstatus.com`
2. **Review GitHub Pages documentation**: `docs.github.com/pages`
3. **Search GitHub Community**: `github.community`
4. **Contact support**: If all else fails, contact GitHub support

## ✅ Deployment Checklist

Before going live:

- [ ] All content reviewed and proofread
- [ ] Images optimized and properly referenced
- [ ] Contact information updated
- [ ] Google AdSense code ready (if using)
- [ ] Custom domain configured (if using)
- [ ] All links tested and working
- [ ] Mobile responsiveness verified
- [ ] Search functionality tested
- [ ] Legal pages completed
- [ ] Analytics/tracking configured
- [ ] Backup of all files created

After deployment:

- [ ] Site loads correctly at GitHub Pages URL
- [ ] Custom domain working (if configured)
- [ ] HTTPS enabled and working
- [ ] All pages accessible
- [ ] Search functionality working
- [ ] Mobile version tested
- [ ] Google Search Console configured
- [ ] Analytics tracking verified

## 🎉 Congratulations!

Your Shell2R website is now live and accessible to the world! You've successfully:

- ✅ Created a professional bioinformatics learning platform
- ✅ Deployed it to GitHub Pages
- ✅ Configured custom domain (if applicable)
- ✅ Set up a system for easy content updates

Your website is now ready to help others learn bioinformatics and single-cell RNA-seq analysis!

---

**Need help with deployment? Have questions about customization? Feel free to reach out for support!**

