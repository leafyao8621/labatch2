#include <gtk/gtk.h>
#include <stdio.h>

static void run(GtkWidget* this, gpointer e) {
    GtkGrid* grid = (GtkGrid*)e;
    const char* ifn = gtk_entry_get_text(GTK_ENTRY(gtk_grid_get_child_at(grid, 2, 0)));
    const char* ofn = gtk_entry_get_text(GTK_ENTRY(gtk_grid_get_child_at(grid, 2, 1)));
    const char* cfn = gtk_entry_get_text(GTK_ENTRY(gtk_grid_get_child_at(grid, 2, 2)));
    FILE* ifs = fopen(ifn, "r");
    FILE* ofs = fopen(ofn, "w");
    FILE* cfs = fopen(cfn, "r");
    char b[80];
    fscanf(ifs, "%s", b);
    fputs("test", ofs);
    fputs(b, ofs);
    char b1[80];
    fscanf(cfs, "%s", b1);
    fputs(b1, ofs);
    fclose(ifs);
    fclose(ofs);
    fclose(cfs);
}

static void activate(GtkApplication* app, gpointer user_data) {
    GtkWidget* window;
    GtkWidget* grid;
    GtkWidget* box;
    GtkWidget* label;
    GtkWidget* entry;
    GtkWidget* label1;
    GtkWidget* entry1;
    GtkWidget* label2;
    GtkWidget* entry2;
    GtkWidget* button;
    window = gtk_application_window_new(app);
    gtk_window_set_title(GTK_WINDOW(window), "lab2");
    gtk_window_set_default_size(GTK_WINDOW(window), 500, 400);
    grid = gtk_grid_new();
    gtk_grid_set_column_spacing(GTK_GRID(grid), 5);
    label = gtk_label_new("Input file name");
    gtk_grid_attach(GTK_GRID(grid), label, 1, 0, 1, 1);
    entry = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), entry, 2, 0, 1, 1);
    label1 = gtk_label_new("Output file name");
    gtk_grid_attach(GTK_GRID(grid), label1, 1, 1, 1, 1);
    entry1 = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), entry1, 2, 1, 1, 1);
    label2 = gtk_label_new("Config file name");
    gtk_grid_attach(GTK_GRID(grid), label2, 1, 2, 1, 1);
    entry2 = gtk_entry_new();
    gtk_grid_attach(GTK_GRID(grid), entry2, 2, 2, 1, 1);
    button = gtk_button_new_with_label("run");
    g_signal_connect(button, "clicked", G_CALLBACK(run), grid);
    gtk_widget_set_halign(grid, GTK_ALIGN_CENTER);
    gtk_widget_set_halign(button, GTK_ALIGN_CENTER);
    box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
    gtk_box_set_homogeneous(GTK_BOX(box), FALSE);
    gtk_box_pack_start(GTK_BOX(box), grid, TRUE, TRUE, 0);
    gtk_box_pack_start(GTK_BOX(box), button, TRUE, TRUE, 0);
    gtk_widget_set_halign(box, GTK_ALIGN_CENTER);
    gtk_widget_set_valign(box, GTK_ALIGN_CENTER);
    gtk_container_add(GTK_CONTAINER(window), box);
    gtk_widget_show_all(window);
}

int main(int argc, char** argv) {
    GtkApplication* app;
    int status;
    app = gtk_application_new("org.gtk.lab2", G_APPLICATION_FLAGS_NONE);
    g_signal_connect(app, "activate", G_CALLBACK(activate), NULL);
    status = g_application_run(G_APPLICATION(app), argc, argv);
    g_object_unref(app);
    return status;
}
