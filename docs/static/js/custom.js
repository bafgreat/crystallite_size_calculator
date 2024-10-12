// custom.js - Enhancements for MkDocs with Material theme

// 1. Scroll to Top Button
document.addEventListener('DOMContentLoaded', function () {
    const scrollToTopBtn = document.createElement('button');
    scrollToTopBtn.id = 'scrollToTopBtn';
    scrollToTopBtn.innerText = '↑';
    scrollToTopBtn.style.cssText = `
        position: fixed;
        bottom: 20px;
        right: 20px;
        display: none;
        background-color: #333;
        color: white;
        border: none;
        border-radius: 5px;
        padding: 10px;
        cursor: pointer;
        z-index: 1000;
    `;

    document.body.appendChild(scrollToTopBtn);

    window.addEventListener('scroll', function () {
        if (window.scrollY > 300) {
            scrollToTopBtn.style.display = 'block';
        } else {
            scrollToTopBtn.style.display = 'none';
        }
    });

    scrollToTopBtn.addEventListener('click', function () {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    });
});

// 2. Smooth Scroll for Anchor Links
document.addEventListener('DOMContentLoaded', function () {
    const links = document.querySelectorAll('a[href^="#"]');

    links.forEach(link => {
        link.addEventListener('click', function (e) {
            e.preventDefault();
            const targetId = link.getAttribute('href');
            const targetElement = document.querySelector(targetId);

            if (targetElement) {
                targetElement.scrollIntoView({ behavior: 'smooth' });
            }
        });
    });
});

// 3. Highlight Active Menu Item in Navigation
document.addEventListener('DOMContentLoaded', function () {
    const currentPath = window.location.pathname;
    const menuItems = document.querySelectorAll('nav a');

    menuItems.forEach(item => {
        if (item.getAttribute('href') === currentPath) {
            item.classList.add('active');
            item.style.fontWeight = 'bold'; // Custom style for active link
        }
    });
});

// 4. Toggle Dark Mode (if supported by theme)
document.addEventListener('DOMContentLoaded', function () {
    const themeSwitcher = document.createElement('button');
    themeSwitcher.id = 'themeSwitcher';
    themeSwitcher.innerText = 'Toggle Dark Mode';
    themeSwitcher.style.cssText = `
        position: fixed;
        bottom: 60px;
        right: 20px;
        background-color: #333;
        color: white;
        border: none;
        border-radius: 5px;
        padding: 10px;
        cursor: pointer;
        z-index: 1000;
    `;

    document.body.appendChild(themeSwitcher);

    themeSwitcher.addEventListener('click', function () {
        const currentTheme = document.documentElement.getAttribute('data-md-color-scheme');
        const newTheme = currentTheme === 'default' ? 'slate' : 'default';
        document.documentElement.setAttribute('data-md-color-scheme', newTheme);
    });
});
