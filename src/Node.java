import java.util.ArrayList;
import java.util.List;

public class Node<T> {
    private List<Node<T>> children = new ArrayList<Node<T>>();
    private Node<T> parent = null;
    private T data = null;
    private String label;
    private String labelValue = "E";
    private int labelIndex = -1;

    public Node(T data) {
        this.data = data;
        this.label = null;
    }
    
    public void setLabelValue(String str){
    	this.labelValue = str;
    }
    
    public String getLabelValue(){
    	return labelValue;
    }
    public void setLabel(String l){
    	this.label = l;
    }
    
    public Node<T> getParent(){
    	return parent;
    }
    
    public String getLabel(){
    	return label;
    }
    public Node(T data, Node<T> parent) {
        this.data = data;
        this.parent = parent;
    }

    public List<Node<T>> getChildren() {
        return children;
    }

    public void setParent(Node<T> parent) {
        this.parent = parent;
    }

    public void addChild(T data) {
        Node<T> child = new Node<T>(data);
        child.setParent(this);
        this.children.add(child);
    }

    public void addChild(Node<T> child) {
        child.setParent(this);
        this.children.add(child);
    }
    
    public void addChild(Node<T> child, String lblValue) {
        child.setParent(this);
        child.setLabelValue(lblValue);
        this.children.add(child);
    }

    public T getData() {
        return this.data;
    }

    public void setData(T data) {
        this.data = data;
    }

    public boolean isRoot() {
        return (this.parent == null);
    }

    public boolean isLeaf() {
        if(this.children.size() == 0) 
            return true;
        else 
            return false;
    }
    
    public String getParentLabel(){
    	return parent.getLabel();
    }

    public void removeParent() {
        this.parent = null;
    }
    
    public boolean hasLabelIndex(){
    	if(labelIndex != -1){
    		return true;
    	}
    	return false;
    }
    
    public void setLabelIndex(int i){
    	this.labelIndex = i;
    }
    
    public int getLabelIndex(){
    	return labelIndex;
    }
}