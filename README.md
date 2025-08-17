# Shell2R - Educational Bioinformatics Website

A professional, responsive website for learning bioinformatics, shell commands, HPC tutorials, and single-cell RNA-seq analysis. Built with pure HTML, CSS, and JavaScript for easy deployment and maintenance.

## 🌟 Features

- **Clean, Modern Design**: Professional layout optimized for educational content
- **Responsive**: Works perfectly on desktop, tablet, and mobile devices
- **Tutorial System**: Automatic tutorial loading and display from Markdown files
- **Search & Filter**: Real-time search and category filtering
- **Code Highlighting**: Beautiful syntax highlighting for R, Python, and Bash
- **Google AdSense Ready**: Pre-integrated ad placement sections
- **SEO Optimized**: Proper meta tags and structured content

## 📁 Project Structure

```
shell2r-simple/
├── index.html              # Main homepage
├── style.css              # Main stylesheet
├── script.js              # JavaScript functionality
├── pages/                 # Static HTML pages
│   ├── about.html         # About page
│   ├── contact.html       # Contact page
│   ├── privacy-policy.html
│   ├── terms-of-service.html
│   ├── disclaimer.html
│   └── cookie-policy.html
├── lessons/               # Tutorial markdown files
│   ├── introduction-to-bioinformatics.md
│   ├── command-line-basics-detailed.md
│   ├── conda-mamba-installation-guide.md
│   └── single-cell-rnaseq-introduction.md
├── images/                # Tutorial images and assets
│   ├── bioinformatics-intro.png
│   ├── command-line-guide.png
│   ├── conda-mamba-guide.png
│   └── single-cell-intro.png
└── docs/                  # Documentation files
    ├── DEPLOYMENT_GUIDE.md
    └── TUTORIAL_GUIDE.md
```

## 🚀 Quick Start

### Local Development

1. **Clone or download** the project files
2. **Open `index.html`** in your web browser to preview the site
3. **For full functionality**, run a local server:
   ```bash
   # Using Python 3
   python -m http.server 8000
   
   # Using Node.js (if you have it installed)
   npx serve .
   
   # Using PHP (if you have it installed)
   php -S localhost:8000
   ```
4. **Open** `http://localhost:8000` in your browser

### GitHub Pages Deployment

1. **Create a new GitHub repository**
2. **Upload all files** to the repository root
3. **Enable GitHub Pages** in repository settings:
   - Go to Settings → Pages
   - Source: "Deploy from a branch"
   - Branch: "main" or "master"
   - Folder: "/ (root)"
4. **Configure custom domain** (optional):
   - Add your domain in the "Custom domain" field
   - Configure DNS A records to point to GitHub Pages

## 📝 Adding New Tutorials

1. **Create a new Markdown file** in the `lessons/` folder:
   ```
   lessons/your-new-tutorial.md
   ```

2. **Add front matter** at the top of your file:
   ```markdown
   ---
   title: "Your Tutorial Title"
   date: "2025-01-17"
   author: "Your Name"
   category: "Bioinformatics"
   excerpt: "Brief description of your tutorial"
   image: "images/your-tutorial-image.png"
   ---
   
   # Your Tutorial Content
   
   Your tutorial content goes here...
   ```

3. **Update `script.js`** to include your new tutorial:
   ```javascript
   const tutorialFiles = [
       'introduction-to-bioinformatics.md',
       'command-line-basics-detailed.md',
       'conda-mamba-installation-guide.md',
       'single-cell-rnaseq-introduction.md',
       'your-new-tutorial.md'  // Add this line
   ];
   ```

4. **Add tutorial image** to the `images/` folder (optional)

5. **Commit and push** to GitHub - your tutorial will appear automatically!

## 🎨 Customization

### Colors and Styling
Edit `style.css` to customize:
- Color scheme (CSS custom properties at the top)
- Typography and fonts
- Layout and spacing
- Component styling

### Navigation
Update navigation links in:
- `index.html` (main navigation)
- `pages/*.html` (page-specific navigation)

### Content
- **Homepage content**: Edit `index.html`
- **Static pages**: Edit files in `pages/` folder
- **Tutorials**: Add/edit Markdown files in `lessons/` folder

## 🔧 Technical Details

### Browser Compatibility
- Modern browsers (Chrome, Firefox, Safari, Edge)
- Mobile browsers (iOS Safari, Chrome Mobile)
- Progressive enhancement for older browsers

### Performance
- Lightweight CSS and JavaScript
- Optimized images
- Minimal external dependencies
- Fast loading times

### SEO Features
- Semantic HTML structure
- Meta tags for social sharing
- Structured data markup
- Mobile-friendly design

## 📊 Google AdSense Integration

The website is pre-configured for Google AdSense:

1. **Replace placeholder ad code** in relevant files
2. **Add your AdSense publisher ID**
3. **Configure ad placements** as needed
4. **Ensure compliance** with all legal pages included

## 🆘 Troubleshooting

### Tutorials Not Loading
- Check that tutorial files are in the `lessons/` folder
- Verify file names match those listed in `script.js`
- Ensure proper front matter format in Markdown files
- Check browser console for error messages

### Local Development Issues
- Use a local server instead of opening files directly
- Check for CORS issues when loading Markdown files
- Ensure all file paths are correct and case-sensitive

### GitHub Pages Issues
- Verify all files are uploaded to repository root
- Check that GitHub Pages is enabled in repository settings
- Allow 5-10 minutes for changes to propagate
- Check repository Actions tab for build errors

## 📞 Support

For questions or issues:
- Check the `TUTORIAL_GUIDE.md` for detailed tutorial creation instructions
- Review `DEPLOYMENT_GUIDE.md` for deployment help
- Open an issue in the GitHub repository

## 📄 License

This project is open source and available under the MIT License.

---

**Shell2R** - Learn Bioinformatics and Single-cell RNA-seq Step by Step

