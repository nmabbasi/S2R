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
        const tutorialFiles = [
            'introduction-to-bioinformatics.md',
            'command-line-basics-detailed.md',
            'conda-mamba-installation-guide.md',
            'single-cell-rnaseq-introduction.md'
        ];
        
        tutorials = [];
        
        for (const file of tutorialFiles) {
            try {
                const response = await fetch(`lessons/${file}`);
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
        
        updateTutorialsList();
        updateSidebar();
        
    } catch (error) {
        console.error('Error loading tutorials:', error);
        showFallbackTutorials();
    }
}

// Parse tutorial markdown content
function parseTutorial(content, filename) {
    try {
        const frontMatterMatch = content.match(/^---\s*\n([\s\S]*?)\n---\s*\n([\s\S]*)$/);
        if (!frontMatterMatch) {
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
        const mainContent = frontMatterMatch[9];
        const metadata = {};
        frontMatter.split('\n').forEach(line => {
            const match = line.match(/^(\w+):\s*(.+)$/);
            if (match) {
                metadata[match[10]] = match[9].replace(/^["']|["']$/g, '');
            }
        });
        let excerpt = metadata.excerpt || '';
        if (!excerpt) {
            const contentWithoutHeaders = mainContent.replace(/^#+\s+.*$/gm, '');
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

// Show fallback tutorials when markdown files can't be loaded (kept for robustness)
function showFallbackTutorials() {
    // Fallback content omitted for brevity...
    // You may retain or remove as needed.
}

// Update tutorials list on homepage
function updateTutorialsList() {
    const tutorialsList = document.getElementById('tutorials-list');
    if (!tutorialsList) return;
    
    if (tutorials.length === 0) {
        tutorialsList.innerHTML = '<p class="text-gray-500">Loading tutorials...</p>';
        return;
    }
    
    const latestTutorials = tutorials.slice(0, 5);
    
    tutorialsList.innerHTML = latestTutorials.map(tutorial => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${tutorial.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${tutorial.category}</span>
                <span class="text-sm text-gray-500">${formatDate(tutorial.date)}</span>
            </div>
            <h3>${tutorial.title}</h3>
            <div class="meta">by ${tutorial.author}</div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">Read More →</a>
        </article>
    `).join('');
}

// Update sidebar with categories and recent posts
function updateSidebar() {
    updateCategories();
    updateRecentPosts();
}

function updateCategories() {
    const categoriesList = document.getElementById('categories-list');
    if (!categoriesList) return;
    
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

    const tutorialContent = document.getElementById('tutorial-content');
    tutorialContent.innerHTML = `
        <article class="tutorial-content bg-white rounded-lg shadow-md p-8">
            <header class="mb-8">
                <div class="flex items-center justify-between mb-4">
                    <span class="category">${tutorial.category}</span>
                    <button onclick="showTutorials()" class="text-blue-600 hover:text-blue-800">← Back to Tutorials</button>
                </div>
                <h1>${tutorial.title}</h1>
                <div class="flex items-center text-gray-600 text-sm">
                    <span>by ${tutorial.author}</span>
                    <span class="mx-2">•</span>
                    <span>${formatDate(tutorial.date)}</span>
                </div>
            </header>
            <div class="prose max-w-none">${marked.parse(tutorial.content)}</div>
        </article>
    `;
    setTimeout(() => { Prism.highlightAll(); }, 100);
    window.scrollTo(0, 0);
}

function showPage(pageName) {
    hideAllPages();
    document.getElementById('static-page').classList.remove('hidden');
    currentPage = pageName;
    loadStaticPage(pageName);
}

// Load static page content from the same folder as index.html
async function loadStaticPage(pageName) {
    const staticContent = document.getElementById('static-content');
    
    try {
        const response = await fetch(`${pageName}.html`);
        if (response.ok) {
            const content = await response.text();
            staticContent.innerHTML = content;
        } else {
            staticContent.innerHTML = '<p>Page not found.</p>';
        }
    } catch (error) {
        console.error('Error loading page:', error);
        staticContent.innerHTML = '<p>Unable to load page content.</p>';
    }
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
    // You may add code here to highlight active nav link based on `activePage`
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

// Filter by category
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
            <div class="meta">by ${tutorial.author}</div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">Read More →</a>
        </article>
    `).join('');
}

// Search functionality
function handleSearch(event) {
    const query = event.target.value.toLowerCase().trim();
    
    if (query === '') {
        if (currentPage === 'tutorials') {
            updateAllTutorialsList();
        } else {
            updateTutorialsList();
        }
        return;
    }
    
    const filteredTutorials = tutorials.filter(tutorial => 
        tutorial.title.toLowerCase().includes(query) ||
        tutorial.content.toLowerCase().includes(query) ||
        tutorial.category.toLowerCase().includes(query) ||
        tutorial.excerpt.toLowerCase().includes(query)
    );
    
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
                    <div class="meta">by ${tutorial.author}</div>
                    <p class="excerpt">${highlightSearchTerm(tutorial.excerpt, query)}</p>
                    <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">Read More →</a>
                </article>
            `).join('');
        }
    }
}

function highlightSearchTerm(text, term) {
    if (!term) return text;
    const regex = new RegExp(`(${term})`, 'gi');
    return text.replace(regex, '<span class="search-highlight">$1</span>');
}

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
            <div class="meta">by ${tutorial.author}</div>
            <p class="excerpt">${tutorial.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation(); showTutorial('${tutorial.id}')">Read More →</a>
        </article>
    `).join('');
}

function updateCategoryFilter() {
    const categoryFilter = document.getElementById('category-filter');
    if (!categoryFilter) return;
    
    const categories = [...new Set(tutorials.map(t => t.category))].sort();
    
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
