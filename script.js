// Shell2R Website JavaScript
let tutorials = [];
let currentPage = 'home';
let currentTutorial = null;

document.addEventListener('DOMContentLoaded', function() {
    loadTutorials();
    setupEventListeners();
    showHome();
    createBackToTopButton();
});

function setupEventListeners() {
    const searchInput = document.getElementById('search-input');
    if (searchInput) {
        searchInput.addEventListener('input', debounce(handleSearch, 300));
    }
    window.addEventListener('scroll', handleScroll);
}

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
                    if (tutorial) tutorials.push(tutorial);
                }
            } catch (e) {
                console.warn(`Could not load tutorial: ${file}`, e);
            }
        }
        tutorials.sort((a, b) => new Date(b.date) - new Date(a.date));
        updateTutorialsList();
        updateSidebar();
    } catch (error) {
        console.error('Error loading tutorials:', error);
        // fallback logic if needed
    }
}

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
        const metaText = frontMatterMatch[1];
        const mainContent = frontMatterMatch[11];
        const metadata = {};
        metaText.split('\n').forEach(line => {
            const match = line.match(/^(\w+):\s*(.+)$/);
            if (match) metadata[match[12]] = match[11].replace(/^["']|["']$/g, '');
        });
        const excerpt = metadata.excerpt || mainContent.replace(/^#+\s+.*$/gm, '').substring(0, 200).trim() + '...';

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

function updateTutorialsList() {
    const tutorialsList = document.getElementById('tutorials-list');
    if (!tutorialsList) return;
    if (tutorials.length === 0) {
        tutorialsList.innerHTML = '<p class="text-gray-500">Loading tutorials...</p>';
        return;
    }
    const latestFive = tutorials.slice(0, 5);
    tutorialsList.innerHTML = latestFive.map(t => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${t.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${t.category}</span>
                <span class="text-sm text-gray-500">${formatDate(t.date)}</span>
            </div>
            <h3>${t.title}</h3>
            <div class="meta">by ${t.author}</div>
            <p class="excerpt">${t.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation();showTutorial('${t.id}')">Read More →</a>
        </article>
    `).join('');
}

function updateSidebar() {
    updateCategories();
    updateRecentPosts();
}

function updateCategories() {
    const categoriesList = document.getElementById('categories-list');
    if (!categoriesList) return;
    const count = {};
    tutorials.forEach(t => count[t.category] = (count[t.category] || 0) + 1);
    categoriesList.innerHTML = Object.keys(count).sort().map(cat => `
        <div class="flex items-center justify-between">
            <a href="#" class="sidebar-link" onclick="filterByCategory('${cat}')">${cat}</a>
            <span class="category-count">${count[cat]}</span>
        </div>
    `).join('');
}

function updateRecentPosts() {
    const recentPosts = document.getElementById('recent-posts');
    if (!recentPosts) return;
    const recent = tutorials.slice(0, 5);
    recentPosts.innerHTML = recent.map(t => `
        <div class="cursor-pointer hover:bg-gray-50 p-2 rounded" onclick="showTutorial('${t.id}')">
            <h4 class="text-sm font-medium text-gray-900 line-clamp-2 mb-1">${t.title}</h4>
            <p class="text-xs text-gray-500">${formatDate(t.date)}</p>
        </div>
    `).join('');
}

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

function showTutorial(id) {
    const tutorial = tutorials.find(t => t.id === id);
    if (!tutorial) return;
    hideAllPages();
    document.getElementById('tutorial-page').classList.remove('hidden');
    currentPage = 'tutorial';
    currentTutorial = tutorial;
    const container = document.getElementById('tutorial-content');
    container.innerHTML = `
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
    setTimeout(() => Prism.highlightAll(), 100);
    window.scrollTo(0, 0);
}

function showPage(pageName) {
    hideAllPages();
    document.getElementById('static-page').classList.remove('hidden');
    currentPage = pageName;
    loadStaticPage(pageName);
}

// Load static page content from files in same folder as index.html
async function loadStaticPage(name) {
    const container = document.getElementById('static-content');
    try {
        const response = await fetch(`${name}.html`);
        if (response.ok) {
            const html = await response.text();
            container.innerHTML = html;
        } else {
            container.innerHTML = '<p>Page not found.</p>';
        }
    } catch (e) {
        console.error('Error loading static page:', e);
        container.innerHTML = '<p>Unable to load content.</p>';
    }
}

function hideAllPages() {
    document.querySelectorAll('.page-content').forEach(p => p.classList.add('hidden'));
}

function updateActiveNavLink(page) {
    document.querySelectorAll('.nav-link').forEach(link => {
        link.classList.remove('text-blue-600');
        link.classList.add('text-gray-700');
    });
    // Optionally highlight active nav link e.g., with data-attributes
}

function formatDate(d) {
    return new Date(d).toLocaleDateString('en-US', { year: 'numeric', month: 'long', day: 'numeric' });
}

function debounce(fn, delay) {
    let timeout;
    return function(...args) {
        clearTimeout(timeout);
        timeout = setTimeout(() => fn.apply(this, args), delay);
    };
}

function toggleMobileMenu() {
    const menu = document.getElementById('mobile-menu');
    menu.classList.toggle('hidden');
}

function createBackToTopButton() {
    const btn = document.createElement('div');
    btn.className = 'back-to-top';
    btn.innerHTML = '↑';
    btn.onclick = () => window.scrollTo({ top: 0, behavior: 'smooth' });
    document.body.appendChild(btn);
}

function handleScroll() {
    const btn = document.querySelector('.back-to-top');
    if (window.pageYOffset > 300) btn.classList.add('visible');
    else btn.classList.remove('visible');
}

function filterByCategory(category) {
    const buttons = document.querySelectorAll('.category-btn');
    buttons.forEach(btn => btn.classList.remove('active'));
    const activeBtn = document.querySelector(`[data-category="${category}"]`);
    if (activeBtn) activeBtn.classList.add('active');
    
    const list = document.getElementById('all-tutorials-list');
    if (!list) return;
    const filtered = category === 'all' ? tutorials : tutorials.filter(t => t.category === category);
    list.innerHTML = filtered.map(t => `
        <article class="tutorial-card cursor-pointer" onclick="showTutorial('${t.id}')">
            <div class="flex items-start justify-between mb-3">
                <span class="category">${t.category}</span>
                <span class="text-sm text-gray-500">${formatDate(t.date)}</span>
            </div>
            <h3>${t.title}</h3>
            <div class="meta">by ${t.author}</div>
            <p class="excerpt">${t.excerpt}</p>
            <a href="#" class="read-more" onclick="event.stopPropagation();showTutorial('${t.id}')">Read More →</a>
        </article>
    `).join('');
}

function handleSearch(e) {
    const query = e.target.value.toLowerCase().trim();
    if (!query) {
        if (currentPage === 'tutorials') updateAllTutorialsList();
        else updateTutorialsList();
        return;
    }
    const filtered = tutorials.filter(t =>
        t.title.toLowerCase().includes(query) ||
        t.content.toLowerCase().includes(query) ||
        t.category.toLowerCase().includes(query) ||
        t.excerpt.toLowerCase().includes(query)
    );
    const targetList = currentPage === 'tutorials' ? document.getElementById('all-tutorials-list') : document.getElementById('tutorials-list');
    if (!targetList) return;
    if (filtered.length === 0) {
        targetList.innerHTML = '<p class="text-gray-500 text-center py-8">No tutorials found.</p>';
    } else {
        targetList.innerHTML = filtered.map(t => `
            <article class="tutorial-card cursor-pointer" onclick="showTutorial('${t.id}')">
                <div class="flex items-start justify-between mb-3">
                    <span class="category">${t.category}</span>
                    <span class="text-sm text-gray-500">${formatDate(t.date)}</span>
                </div>
                <h3>${highlightSearchTerm(t.title, query)}</h3>
                <div class="meta">by ${t.author}</div>
                <p class="excerpt">${highlightSearchTerm(t.excerpt, query)}</p>
                <a href="#" class="read-more" onclick="event.stopPropagation();showTutorial('${t.id}')">Read More →</a>
            </article>
        `).join('');
    }
}

function highlightSearchTerm(text, term) {
    if (!term) return text;
    const regex = new RegExp(`(${term})`, 'gi');
    return text.replace(regex, '<span class="search-highlight">$1</span>');
}

// Implement other helper functions as needed...

