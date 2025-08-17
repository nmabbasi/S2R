// Shell2R Website JavaScript

// Global variables
let tutorials = [];
let currentPage = 'home';
let currentTutorial = null;

// Initialize the website
document.addEventListener('DOMContentLoaded', function() {
    loadTutorials();
    setupEventListeners();
    showHome();
    createBackToTopButton();
});

// Setup event listeners
function setupEventListeners() {
    // Search functionality
    const searchInput = document.getElementById('search-input');
    if (searchInput) {
        searchInput.addEventListener('input', debounce(handleSearch, 300));
    }
    
    // Back to top button
    window.addEventListener('scroll', handleScroll);
}

// Load tutorials from markdown files
async function loadTutorials() {
    try {
        // List of tutorial files (you'll add more as you create them)
        const tutorialFiles = [
            'introduction-to-bioinformatics.md',
            'command-line-basics-detailed.md',
            'conda-mamba-installation-guide.md',
            'single-cell-rnaseq-introduction.md'
        ];
        
        tutorials = [];
        
        for (const file of tutorialFiles) {
            try {
                const response = await fetch(`./lessons/${file}`);
                if (response.ok) {
                    const content = await response.text();
                    const tutorial = parseTutorial(content, file);
                    if (tutorial) {
                        tutorials.push(tutorial);
                    }
                }
            } catch (error) {
                console.warn(`Could not load tutorial: ${file}`, error);
            }
        }
        
        // Sort tutorials by date (newest first)
        tutorials.sort((a, b) => new Date(b.date) - new Date(a.date));
        
        // Update the UI
        updateTutorialsList();
        updateSidebar();
        
    } catch (error) {
        console.error('Error loading tutorials:', error);
        // Show fallback content
        showFallbackTutorials();
    }
}

// Parse tutorial markdown content
function parseTutorial(content, filename) {
    try {
        // Extract front matter (metadata between --- lines)
        const frontMatterMatch = content.match(/^---\s*\n([\s\S]*?)\n---\s*\n([\s\S]*)$/);
        
        if (!frontMatterMatch) {
            // If no front matter, create basic tutorial info
            return {
                id: filename.replace('.md', ''),
                title: filename.replace('.md', '').replace(/-/g, ' ').replace(/\b\w/g, l => l.toUpperCase()),
                date: new Date().toISOString().split('T')[0],
                author: 'Shell2R Team',
                category: 'Bioinformatics',
                excerpt: content.substring(0, 200) + '...',
                content: content,
                filename: filename
            };
        }
        
        const frontMatter = frontMatterMatch[1];
        const mainContent = frontMatterMatch[2];
        
        // Parse front matter
        const metadata = {};
        frontMatter.split('\n').forEach(line => {
            const match = line.match(/^(\w+):\s*(.+)$/);
            if (match) {
                metadata[match[1]] = match[2].replace(/^["']|["']$/g, '');
            }
        });
        
        // Extract excerpt from content if not provided
        let excerpt = metadata.excerpt || '';
        if (!excerpt) {
            const contentWithoutHeaders = mainContent.replace(/^#+\s+.*/gm, '');
            excerpt = contentWithoutHeaders.substring(0, 200).trim() + '...';
        }
        
        return {
            id: filename.replace('.md', ''),
            title: metadata.title || 'Untitled Tutorial',
            date: metadata.date || new Date().toISOString().split('T')[0],
            author: metadata.author || 'Shell2R Team',
            category: metadata.category || 'Bioinformatics',
            excerpt: excerpt,
            content: mainContent,
            filename: filename
        };
        
    } catch (error) {
        console.error('Error parsing tutorial:', filename, error);
        return null;
    }
}

// Show fallback tutorials when markdown files can't be loaded
function showFallbackTutorials() {
    tutorials = [
        {
            id: 'introduction-to-bioinformatics',
            title: 'Introduction to Bioinformatics: Getting Started with Biological Data Analysis',
            date: '2024-08-15',
            author: 'Shell2R Team',
            category: 'Bioinformatics',
            excerpt: 'Learn the fundamentals of bioinformatics and discover how computational methods are revolutionizing biological research. This tutorial covers basic concepts, tools, and workflows.',
            content: `# Introduction to Bioinformatics

## What is Bioinformatics?

Bioinformatics is an interdisciplinary field that combines biology, computer science, mathematics, and statistics to analyze and interpret biological data. With the explosion of biological data from genomics, proteomics, and other high-throughput technologies, bioinformatics has become essential for modern biological research.

## Why Learn Bioinformatics?

In today's data-driven world, biological research generates massive amounts of information. Consider these statistics:

- The human genome contains approximately 3.2 billion base pairs
- A single RNA-seq experiment can generate millions of sequencing reads
- Protein databases contain information on hundreds of thousands of proteins

Without computational tools, analyzing this data would be impossible. Bioinformatics enables researchers to:

- Process large datasets efficiently
- Identify patterns in biological data
- Make predictions about biological functions
- Accelerate discovery in medicine and biology

## Getting Started

To begin your bioinformatics journey, you'll need to master several key areas:

1. **Command Line Skills** - Essential for running bioinformatics tools
2. **Programming** - R and Python are the most popular languages
3. **Statistics** - Understanding data analysis and interpretation
4. **Biology** - Domain knowledge is crucial for meaningful analysis

## Next Steps

Ready to dive deeper? Check out our other tutorials on command line basics, package management with Conda, and single-cell RNA-seq analysis.`,
            filename: 'introduction-to-bioinformatics.md'
        },
        {
            id: 'command-line-basics',
            title: 'Command Line Mastery: A Detailed Guide for Bioinformatics Beginners',
            date: '2024-08-14',
            author: 'Shell2R Team',
            category: 'Shell Commands',
            excerpt: 'Master the command line from scratch — learn essential Unix commands, file manipulation, and text processing skills that every bioinformatician needs to succeed.',
            content: `# Command Line Mastery for Bioinformatics

## Why the Command Line Matters

The command line is like learning to drive a manual transmission car. Sure, automatic is easier to start with, but once you master manual, you have complete control over the machine. In bioinformatics, that control translates to:

- **Processing massive datasets** that would crash graphical programs
- **Automating repetitive tasks** that would take hours manually
- **Connecting tools together** in powerful workflows
- **Working on remote servers** where GUIs aren't available

## Essential Commands

### Navigation
\`\`\`bash
pwd          # Print working directory
ls           # List files
ls -la       # List all files with details
cd           # Change directory
cd ..        # Go up one level
cd ~         # Go to home directory
\`\`\`

### File Operations
\`\`\`bash
cp file1 file2       # Copy file
mv file1 file2       # Move/rename file
rm file              # Remove file
mkdir directory      # Create directory
rmdir directory      # Remove empty directory
\`\`\`

### Text Processing
\`\`\`bash
cat file.txt         # Display file content
head file.txt        # Show first 10 lines
tail file.txt        # Show last 10 lines
grep "pattern" file  # Search for pattern
wc -l file.txt       # Count lines
\`\`\`

## Bioinformatics Examples

### Count sequences in a FASTA file
\`\`\`bash
grep -c ">" sequences.fasta
\`\`\`

### Extract sequence IDs
\`\`\`bash
grep ">" sequences.fasta | sed 's/>//'
\`\`\`

### Calculate sequence lengths
\`\`\`bash
awk '/^>/ {if (seq) print length(seq); seq=""; next} {seq=seq$0} END {print length(seq)}' sequences.fasta
\`\`\`

## Conclusion

The command line is your gateway to powerful bioinformatics analysis. Practice these commands regularly, and you'll soon find yourself working more efficiently than ever before.`,
            filename: 'command-line-basics.md'
        },
        {
            id: 'conda-mamba-guide',
            title: 'Conda and Mamba: The Complete Installation and Usage Guide for Bioinformatics',
            date: '2024-08-13',
            author: 'Shell2R Team',
            category: 'Conda',
            excerpt: 'Master package management in bioinformatics with Conda and Mamba — learn installation, environment management, and how to install essential tools like Seurat for single-cell analysis.',
            content: `# Conda and Mamba for Bioinformatics

## Why Package Management Matters

If you've ever spent hours trying to install a bioinformatics tool only to run into dependency conflicts, version mismatches, or the dreaded "it works on my machine" problem — you're not alone. Package management is one of the biggest pain points for researchers entering computational biology.

That's where Conda and Mamba come in. Think of them as your personal assistants for managing software installations — they handle all the messy details of dependencies, versions, and compatibility so you can focus on your research.

## What Are Conda and Mamba?

**Conda** is a package manager and environment management system that was originally created for Python but has evolved to support packages from any language. It's like having a smart librarian who not only knows where every book is but also ensures that when you check out a book, all the related materials you need are available and compatible.

**Mamba** is a reimplementation of Conda that's significantly faster — we're talking about going from minutes to seconds for complex installations. It's essentially Conda with a turbo engine.

## Installing Conda

### Option 1: Miniconda (Recommended)

\`\`\`bash
# Download Miniconda for Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Make it executable
chmod +x Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh
\`\`\`

## Installing Mamba

\`\`\`bash
conda install -c conda-forge mamba
\`\`\`

## Creating Environments

\`\`\`bash
# Create environment for single-cell analysis
mamba create -n single-cell python=3.9

# Activate the environment
conda activate single-cell

# Install Seurat and dependencies
mamba install -c conda-forge -c bioconda r-seurat r-ggplot2 r-dplyr
\`\`\`

## Best Practices

1. **One Environment Per Project**
2. **Document Your Environments**
3. **Pin Important Versions**
4. **Regular Maintenance**

With proper package management, you'll never have to worry about "dependency hell" again!`,
            filename: 'conda-mamba-guide.md'
        },
        {
            id: 'single-cell-intro',
            title: 'Single-cell RNA-seq Analysis: From Raw Data to Biological Insights',
            date: '2024-08-12',
            author: 'Shell2R Team',
            category: 'Single-cell RNA-seq',
            excerpt: 'Discover the revolutionary world of single-cell RNA sequencing — learn how this technology is transforming our understanding of cellular heterogeneity and development.',
            content: `# Single-cell RNA-seq Analysis

## Introduction to Single-cell RNA-seq

Single-cell RNA sequencing (scRNA-seq) is a revolutionary technology that allows us to measure gene expression in individual cells rather than bulk tissue samples. This approach has transformed our understanding of cellular heterogeneity, development, and disease.

## Why Single-cell?

Traditional bulk RNA-seq provides an average expression profile across all cells in a sample, potentially masking important biological differences between cell types or states. Single-cell RNA-seq overcomes this limitation by:

- Revealing cellular heterogeneity within tissues
- Identifying rare cell types and subtypes
- Tracking developmental trajectories
- Understanding cell state transitions
- Discovering new biological mechanisms

## Key Concepts

### Cell Types vs. Cell States
- **Cell types**: Distinct cellular identities (e.g., neurons, T cells, fibroblasts)
- **Cell states**: Temporary conditions within a cell type (e.g., activated, resting, stressed)

### Technical vs. Biological Variation
- **Technical variation**: Introduced by the experimental protocol
- **Biological variation**: True differences between cells

## Analysis Workflow

### 1. Quality Control
\`\`\`r
# Load libraries
library(Seurat)
library(ggplot2)

# Load data
data <- Read10X(data.dir = "filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = data, project = "scRNA_analysis")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
\`\`\`

### 2. Normalization and Scaling
\`\`\`r
# Normalize data
seurat_obj <- NormalizeData(seurat_obj)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale data
seurat_obj <- ScaleData(seurat_obj)
\`\`\`

### 3. Dimensionality Reduction
\`\`\`r
# Principal Component Analysis
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
\`\`\`

### 4. Clustering
\`\`\`r
# Find neighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)

# Find clusters
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Visualize clusters
DimPlot(seurat_obj, reduction = "umap")
\`\`\`

## Common Challenges

1. **Dropout Events**: Not all genes are detected in every cell
2. **Batch Effects**: Technical variation between experiments
3. **Cell Cycle Effects**: Cells in different phases of division
4. **Doublets**: Two cells captured together

## Best Practices

- Always perform thorough quality control
- Use appropriate normalization methods
- Validate findings with independent datasets
- Consider biological context in interpretation

## Conclusion

Single-cell RNA-seq is a powerful technology that continues to evolve rapidly. By understanding the key concepts and following best practices, you can unlock valuable biological insights from your data.

## Next Steps

- Practice with public datasets
- Learn advanced analysis techniques
- Explore trajectory inference methods
- Study cell-cell communication analysis`,
            filename: 'single-cell-intro.md'
        }
    ];
    
    updateTutorialsList();
    updateSidebar();
}

// Update tutorials list on homepage
function updateTutorialsList() {
    const tutorialsList = document.getElementById('tutorials-list');
    if (!tutorialsList) return;
    
    if (tutorials.length === 0) {
        tutorialsList.innerHTML = '<p class="text-gray-500">Loading tutorials...</p>';
        return;
    }
    
    const latestTutorials = tutorials.slice(0, 5); // Show latest 5 tutorials
    
    tutorialsList.innerHTML = latestTutorials.map(tutorial => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${tutorial.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${tutorial.category}</span>
                <span class="text-sm text-gray-500">${formatDate(tutorial.date)}</span>
            </div>
            <h3>${tutorial.title}</h3>
            <div class="meta">
                by ${tutorial.author}
            </div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">
                Read More →
            </a>
        </article>
    `).join('');
}

// Update sidebar with categories and recent posts
function updateSidebar() {
    updateCategories();
    updateRecentPosts();
}

// Update categories list
function updateCategories() {
    const categoriesList = document.getElementById('categories-list');
    if (!categoriesList) return;
    
    // Count tutorials by category
    const categoryCount = {};
    tutorials.forEach(tutorial => {
        categoryCount[tutorial.category] = (categoryCount[tutorial.category] || 0) + 1;
    });
    
    const categories = Object.keys(categoryCount).sort();
    
    categoriesList.innerHTML = categories.map(category => `
        <div class="flex items-center justify-between">
            <a href="#" class="sidebar-link" onclick="filterByCategory('${category}')">${category}</a>
            <span class="category-count">${categoryCount[category]}</span>
        </div>
    `).join('');
}

// Update recent posts
function updateRecentPosts() {
    const recentPosts = document.getElementById('recent-posts');
    if (!recentPosts) return;
    
    const recent = tutorials.slice(0, 5);
    
    recentPosts.innerHTML = recent.map(tutorial => `
        <div class="cursor-pointer hover:bg-gray-50 p-2 rounded" onclick="showTutorial('${tutorial.id}')">
            <h4 class="text-sm font-medium text-gray-900 line-clamp-2 mb-1">${tutorial.title}</h4>
            <p class="text-xs text-gray-500">${formatDate(tutorial.date)}</p>
        </div>
    `).join('');
}

// Navigation functions
function showHome() {
    hideAllPages();
    document.getElementById('home-page').classList.remove('hidden');
    currentPage = 'home';
    updateActiveNavLink('home');
}

function showTutorials() {
    hideAllPages();
    document.getElementById('tutorials-page').classList.remove('hidden');
    currentPage = 'tutorials';
    updateActiveNavLink('tutorials');
    updateAllTutorialsList();
    updateCategoryFilter();
}

function showTutorial(tutorialId) {
    const tutorial = tutorials.find(t => t.id === tutorialId);
    if (!tutorial) return;
    
    hideAllPages();
    document.getElementById('tutorial-page').classList.remove('hidden');
    currentPage = 'tutorial';
    currentTutorial = tutorial;
    
    // Render tutorial content
    const tutorialContent = document.getElementById('tutorial-content');
    tutorialContent.innerHTML = `
        <article class="tutorial-content bg-white rounded-lg shadow-md p-8">
            <header class="mb-8">
                <div class="flex items-center justify-between mb-4">
                    <span class="category">${tutorial.category}</span>
                    <button onclick="showTutorials()" class="text-blue-600 hover:text-blue-800">
                        ← Back to Tutorials
                    </button>
                </div>
                <h1>${tutorial.title}</h1>
                <div class="flex items-center text-gray-600 text-sm">
                    <span>by ${tutorial.author}</span>
                    <span class="mx-2">•</span>
                    <span>${formatDate(tutorial.date)}</span>
                </div>
            </header>
            <div class="prose max-w-none">
                ${marked.parse(tutorial.content)}
            </div>
        </article>
    `;
    
    // Highlight code blocks
    setTimeout(() => {
        Prism.highlightAll();
    }, 100);
    
    // Scroll to top
    window.scrollTo(0, 0);
}

function showPage(pageName) {
    hideAllPages();
    document.getElementById('static-page').classList.remove('hidden');
    currentPage = pageName;
    
    loadStaticPage(pageName);
}

// Load static page content
async function loadStaticPage(pageName) {
    const staticContent = document.getElementById('static-content');
    
    try {
        const response = await fetch(`pages/${pageName}.html`);
        if (response.ok) {
            const content = await response.text();
            staticContent.innerHTML = content;
        } else {
            // Fallback content
            staticContent.innerHTML = getFallbackPageContent(pageName);
        }
    } catch (error) {
        console.error('Error loading page:', error);
        staticContent.innerHTML = getFallbackPageContent(pageName);
    }
}

// Get fallback content for static pages
function getFallbackPageContent(pageName) {
    const fallbackContent = {
        'about': `
            <h1>About Shell2R</h1>
            <p>Shell2R is an educational platform dedicated to teaching bioinformatics and computational biology through hands-on tutorials and practical examples.</p>
            
            <h2>Our Mission</h2>
            <p>We believe that learning bioinformatics should be accessible, practical, and engaging. Our tutorials are designed to take you from beginner to advanced practitioner, covering everything from basic command line skills to complex single-cell RNA-seq analysis.</p>
            
            <h2>What You'll Learn</h2>
            <ul>
                <li>Essential command line skills for bioinformatics</li>
                <li>Package management with Conda and Mamba</li>
                <li>Single-cell RNA-seq analysis techniques</li>
                <li>Data visualization and interpretation</li>
                <li>Best practices for reproducible research</li>
            </ul>
            
            <h2>About the Team</h2>
            <p>Shell2R is maintained by a team of bioinformatics researchers and educators passionate about making computational biology accessible to everyone.</p>
        `,
        'contact': `
            <h1>Contact Us</h1>
            <p>We'd love to hear from you! Whether you have questions about our tutorials, suggestions for new content, or just want to say hello, don't hesitate to reach out.</p>
            
            <h2>Get in Touch</h2>
            <p>Email us at: <a href="mailto:contact@shell2r.com" class="text-blue-600 hover:text-blue-800">contact@shell2r.com</a></p>
            
            <h2>Feedback and Suggestions</h2>
            <p>Your feedback helps us improve our content and create tutorials that are truly useful. If you have:</p>
            <ul>
                <li>Suggestions for new tutorial topics</li>
                <li>Questions about existing content</li>
                <li>Bug reports or technical issues</li>
                <li>General feedback about the site</li>
            </ul>
            <p>Please don't hesitate to contact us!</p>
            
            <h2>Contributing</h2>
            <p>Interested in contributing to Shell2R? We welcome contributions from the community. Contact us to learn more about how you can help make bioinformatics education better for everyone.</p>
        `,
        'privacy-policy': `
            <h1>Privacy Policy</h1>
            <p><strong>Last updated:</strong> August 15, 2024</p>
            
            <h2>Information We Collect</h2>
            <p>We collect information you provide directly to us, such as when you contact us via email.</p>
            
            <h2>How We Use Your Information</h2>
            <p>We use the information we collect to:</p>
            <ul>
                <li>Respond to your inquiries and provide customer support</li>
                <li>Improve our website and services</li>
                <li>Comply with legal obligations</li>
            </ul>
            
            <h2>Google Analytics</h2>
            <p>We use Google Analytics to understand how visitors interact with our website. This service collects information such as how often users visit the site, what pages they visit, and what other sites they used prior to coming to this site.</p>
            
            <h2>Google AdSense</h2>
            <p>We use Google AdSense to display advertisements. Google may use cookies to serve ads based on your prior visits to our website or other websites.</p>
            
            <h2>Contact Us</h2>
            <p>If you have any questions about this Privacy Policy, please contact us at contact@shell2r.com.</p>
        `,
        'terms-of-service': `
            <h1>Terms of Service</h1>
            <p><strong>Last updated:</strong> August 15, 2024</p>
            
            <h2>Acceptance of Terms</h2>
            <p>By accessing and using Shell2R, you accept and agree to be bound by the terms and provision of this agreement.</p>
            
            <h2>Educational Use</h2>
            <p>This website is intended for educational purposes. The content provided is for informational purposes only and should not be considered as professional advice.</p>
            
            <h2>Content License</h2>
            <p>The educational content on this website is provided under a Creative Commons license. You are free to use, share, and adapt the content for educational purposes with proper attribution.</p>
            
            <h2>Disclaimer</h2>
            <p>The information on this website is provided on an "as is" basis. We make no warranties, expressed or implied, and hereby disclaim all other warranties.</p>
            
            <h2>Contact</h2>
            <p>If you have any questions about these Terms of Service, please contact us at contact@shell2r.com.</p>
        `,
        'disclaimer': `
            <h1>Disclaimer</h1>
            <p><strong>Last updated:</strong> August 15, 2024</p>
            
            <h2>Educational Purpose</h2>
            <p>Shell2R is an educational platform designed to provide tutorials and resources related to bioinformatics and computational biology. All content is provided for educational and informational purposes only.</p>
            
            <h2>No Professional Advice</h2>
            <p>The information provided on Shell2R does not constitute professional, scientific, or technical advice. Always consult with qualified professionals for specific research needs.</p>
            
            <h2>Accuracy and Completeness</h2>
            <p>While we strive to provide accurate and up-to-date information, we make no representations or warranties about the completeness, accuracy, or reliability of any information.</p>
            
            <h2>Limitation of Liability</h2>
            <p>In no event shall Shell2R be liable for any direct, indirect, incidental, or consequential damages arising from the use of this website or its content.</p>
            
            <h2>User Responsibility</h2>
            <p>Users are solely responsible for verifying the accuracy and applicability of information before using it in their research or work.</p>
        `,
        'cookie-policy': `
            <h1>Cookie Policy</h1>
            <p><strong>Last updated:</strong> August 15, 2024</p>
            
            <h2>What Are Cookies?</h2>
            <p>Cookies are small text files that are stored on your computer or mobile device when you visit a website. They help websites work more efficiently and provide information to website owners.</p>
            
            <h2>How We Use Cookies</h2>
            <p>We use cookies to:</p>
            <ul>
                <li>Remember your preferences and settings</li>
                <li>Understand how you use our website</li>
                <li>Improve your browsing experience</li>
                <li>Display relevant advertisements</li>
            </ul>
            
            <h2>Types of Cookies We Use</h2>
            <h3>Essential Cookies</h3>
            <p>These cookies are necessary for the website to function properly.</p>
            
            <h3>Analytics Cookies</h3>
            <p>We use Google Analytics to understand how visitors interact with our website.</p>
            
            <h3>Advertising Cookies</h3>
            <p>We use Google AdSense to display relevant advertisements.</p>
            
            <h2>Managing Cookies</h2>
            <p>You can control cookies through your browser settings. However, disabling cookies may affect the functionality of our website.</p>
            
            <h2>Contact Us</h2>
            <p>If you have questions about our use of cookies, please contact us at contact@shell2r.com.</p>
        `
    };
    
    return fallbackContent[pageName] || '<h1>Page Not Found</h1><p>The requested page could not be found.</p>';
}

// Update all tutorials list (for tutorials page)
function updateAllTutorialsList() {
    const allTutorialsList = document.getElementById('all-tutorials-list');
    if (!allTutorialsList) return;
    
    allTutorialsList.innerHTML = tutorials.map(tutorial => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${tutorial.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${tutorial.category}</span>
                <span class="text-sm text-gray-500">${formatDate(tutorial.date)}</span>
            </div>
            <h3>${tutorial.title}</h3>
            <div class="meta">
                by ${tutorial.author}
            </div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">
                Read More →
            </a>
        </article>
    `).join('');
}

// Update category filter buttons
function updateCategoryFilter() {
    const categoryFilter = document.getElementById('category-filter');
    if (!categoryFilter) return;
    
    const categories = [...new Set(tutorials.map(t => t.category))].sort();
    
    // Keep the "All" button and add category buttons
    const allButton = categoryFilter.querySelector('[data-category="all"]');
    categoryFilter.innerHTML = '';
    categoryFilter.appendChild(allButton);
    
    categories.forEach(category => {
        const button = document.createElement('button');
        button.className = 'category-btn bg-white text-gray-700 px-4 py-2 rounded-lg border';
        button.setAttribute('data-category', category);
        button.textContent = category;
        button.onclick = () => filterByCategory(category);
        categoryFilter.appendChild(button);
    });
}

// Filter tutorials by category
function filterByCategory(category) {
    const buttons = document.querySelectorAll('.category-btn');
    buttons.forEach(btn => btn.classList.remove('active'));
    
    const activeButton = document.querySelector(`[data-category="${category}"]`);
    if (activeButton) {
        activeButton.classList.add('active');
    }
    
    const allTutorialsList = document.getElementById('all-tutorials-list');
    if (!allTutorialsList) return;
    
    const filteredTutorials = category === 'all' 
        ? tutorials 
        : tutorials.filter(t => t.category === category);
    
    allTutorialsList.innerHTML = filteredTutorials.map(tutorial => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${tutorial.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${tutorial.category}</span>
                <span class="text-sm text-gray-500">${formatDate(tutorial.date)}</span>
            </div>
            <h3>${tutorial.title}</h3>
            <div class="meta">
                by ${tutorial.author}
            </div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">
                Read More →
            </a>
        </article>
    `).join('');
}

// Search functionality
function handleSearch(event) {
    const query = event.target.value.toLowerCase().trim();
    
    if (query === '') {
        // Show all tutorials
        if (currentPage === 'tutorials') {
            updateAllTutorialsList();
        } else {
            updateTutorialsList();
        }
        return;
    }
    
    // Filter tutorials based on search query
    const filteredTutorials = tutorials.filter(tutorial => 
        tutorial.title.toLowerCase().includes(query) ||
        tutorial.content.toLowerCase().includes(query) ||
        tutorial.category.toLowerCase().includes(query) ||
        tutorial.excerpt.toLowerCase().includes(query)
    );
    
    // Update the appropriate list
    const targetList = currentPage === 'tutorials' 
        ? document.getElementById('all-tutorials-list')
        : document.getElementById('tutorials-list');
    
    if (targetList) {
        if (filteredTutorials.length === 0) {
            targetList.innerHTML = '<p class="text-gray-500 text-center py-8">No tutorials found matching your search.</p>';
        } else {
            targetList.innerHTML = filteredTutorials.map(tutorial => `
                <article class="tutorial-card cursor-pointer" onclick="showTutorial('${tutorial.id}')">
                    <div class="flex items-start justify-between mb-3">
                        <span class="category">${tutorial.category}</span>
                        <span class="text-sm text-gray-500">${formatDate(tutorial.date)}</span>
                    </div>
                    <h3>${highlightSearchTerm(tutorial.title, query)}</h3>
                    <div class="meta">
                        by ${tutorial.author}
                    </div>
                    <p class="excerpt">${highlightSearchTerm(tutorial.excerpt, query)}</p>
                    <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">
                        Read More →
                    </a>
                </article>
            `).join('');
        }
    }
}

// Highlight search terms
function highlightSearchTerm(text, term) {
    if (!term) return text;
    const regex = new RegExp(`(${term})`, 'gi');
    return text.replace(regex, '<span class="search-highlight">$1</span>');
}

// Utility functions
function hideAllPages() {
    document.querySelectorAll('.page-content').forEach(page => {
        page.classList.add('hidden');
    });
}

function updateActiveNavLink(activePage) {
    document.querySelectorAll('.nav-link').forEach(link => {
        link.classList.remove('text-blue-600');
        link.classList.add('text-gray-700');
    });
    
    // This is a simplified approach - in a real implementation,
    // you might want to add data attributes to identify nav links
}

function formatDate(dateString) {
    const date = new Date(dateString);
    return date.toLocaleDateString('en-US', {
        year: 'numeric',
        month: 'long',
        day: 'numeric'
    });
}

function debounce(func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

// Mobile menu toggle
function toggleMobileMenu() {
    const mobileMenu = document.getElementById('mobile-menu');
    mobileMenu.classList.toggle('hidden');
}

// Back to top functionality
function createBackToTopButton() {
    const backToTop = document.createElement('div');
    backToTop.className = 'back-to-top';
    backToTop.innerHTML = '↑';
    backToTop.onclick = () => window.scrollTo({ top: 0, behavior: 'smooth' });
    document.body.appendChild(backToTop);
}

function handleScroll() {
    const backToTop = document.querySelector('.back-to-top');
    if (window.pageYOffset > 300) {
        backToTop.classList.add('visible');
    } else {
        backToTop.classList.remove('visible');
    }
}
